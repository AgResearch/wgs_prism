# 2024 Benjamin J Perry
# MIT License
# Copyright (c) 2024 Benjamin J Perry
# Version: 1.0
# Maintainer: Benjamin J Perry
# Email: ben.perry@agresearch.co.nz


configfile: "config/pipeline_config.yaml"


import os
import pandas as pd

# config dictionary values to be defined on running snakemake with --config flag
for _required in ("RUN", "IN_ROOT", "OUT_ROOT"):
    # Empty is refused as well as missing: pipeline_config.yaml ships these as
    # "", so `in config` would let the default through. A real-looking default
    # would let a forgotten --config silently read, or write over, another run.
    if not config.get(_required):
        raise WorkflowError(f"{_required} must be set, e.g. --config {_required}=...")


rule targets:
    input:
        os.path.join(
            config["OUT_ROOT"],
            "SampleSheet",
            "multiqc",
            f"{config['RUN']}.multiqc.html",
        ),


# Lives outside the bclconvert output dir so it survives Snakemake deleting
# that directory() output on job failure. Note bclconvert_out itself is only
# created as a side effect of preparing the nested reports/fastq_complete
# outputs below (Snakemake pre-creates a directory() output's parent, not the
# directory itself), so the dir already exists (empty) by the time bcl-convert
# runs. Rather than passing --force to tolerate that, the shell below rmdirs
# the empty shells, which turns bcl-convert's own refusal to overwrite into
# the re-demultiplex guard.
bclconvert_logs_persist = os.path.join(config["OUT_ROOT"], "logs", "bclconvert_Logs")


checkpoint run_bclconvert:
    input:
        run_in=os.path.join(config["IN_ROOT"], config["RUN"]),
        sample_sheet=os.path.join(config["OUT_ROOT"], "SampleSheet.csv"),
    output:
        bclconvert_out=directory(
            os.path.join(config["OUT_ROOT"], "SampleSheet", "bclconvert")
        ),
        reports=directory(
            os.path.join(config["OUT_ROOT"], "SampleSheet", "bclconvert", "Reports")
        ),
        fastq_complete=os.path.join(
            config["OUT_ROOT"],
            "SampleSheet",
            "bclconvert",
            "Logs",
            "FastqComplete.txt",
        ),
    log:
        os.path.join(config["OUT_ROOT"], "logs", "run_bclconvert.log"),
    benchmark:
        os.path.join(config["OUT_ROOT"], "benchmarks", "run_bclconvert.txt")
    threads: 24
    # No retry is possible: the guard below refuses to write over existing
    # output, so every restart fails the same way. Without this, the profile's
    # restart-times spends all 5 attempts on a demux that cannot succeed, and
    # attempt-scaled resources never get a chance to apply. One attempt only,
    # sized for the worst case, then a clear message telling the operator what
    # to remove.
    retries: 0
    resources:
        mem_gb=96,
        time=600,
        partition="compute,hugemem,vgpu",
    shell:
        """
        export PATH=/agr/persist/apps/src/b/BCL-Convert:$PATH

        echo "bcl-convert version in use:" > {log}
        bcl-convert -V >> {log} 2>&1
        echo >> {log}

        # bcl-convert refuses to write into a destination that already exists, and we
        # rely on that so a previous demultiplex is never silently overwritten.
        # Snakemake creates bclconvert/ and bclconvert/Logs/ from the declared outputs
        # before this runs, so clear those empty shells first. If real output is
        # present the rmdir fails and the directory survives the check below.
        #
        # This rmdir and the deliberate absence of --force below are a matched
        # pair: drop either one alone and bcl-convert fails on every fresh run,
        # not just on a re-demultiplex.
        rmdir {output.bclconvert_out}/Logs {output.bclconvert_out} 2>/dev/null || true

        if [ -d {output.bclconvert_out} ]; then
            echo "error: {output.bclconvert_out} already contains demultiplexed output; refusing to overwrite. Remove it if you intend to re-demultiplex." >&2
            exit 1
        fi

        # Capture exit status ourselves (Snakemake's `set -e` would otherwise
        # abort before the log-rescue step below runs on failure).
        bcl_status=0
        bcl-convert \
            --bcl-input-directory {input.run_in} \
            --sample-sheet {input.sample_sheet} \
            --output-directory {output.bclconvert_out} >> {log} 2>&1 || bcl_status=$?

        # Rescue bcl-convert's Logs/ before Snakemake deletes the output dir on
        # failure; `ls` guard avoids tripping `set -e` if Logs/ never appeared.
        if ls {output.bclconvert_out}/Logs/*.log >/dev/null 2>&1; then
            mkdir -p {bclconvert_logs_persist}
            cp -f {output.bclconvert_out}/Logs/*.log {bclconvert_logs_persist}/

            echo >> {log}
            echo "===== bcl-convert Logs/ (preserved in {bclconvert_logs_persist}) =====" >> {log}
            cat {output.bclconvert_out}/Logs/*.log >> {log}
        fi

        exit $bcl_status
        """


fastqc_out_root = os.path.join(
    config["OUT_ROOT"], "SampleSheet", "fastqc_run", "fastqc"
)


rule run_fastqc:
    input:
        fastq=os.path.join(
            config["OUT_ROOT"], "SampleSheet", "bclconvert", "{sample}.fastq.gz"
        ),
    output:
        zip=os.path.join(
            config["OUT_ROOT"],
            "SampleSheet",
            "fastqc_run",
            "fastqc",
            "{sample}_fastqc.zip",
        ),
        html=os.path.join(
            config["OUT_ROOT"],
            "SampleSheet",
            "fastqc_run",
            "fastqc",
            "{sample}_fastqc.html",
        ),
    log:
        os.path.join(config["OUT_ROOT"], "logs", "run_fastqc.{sample}.log"),
    conda:
        "envs/fastqc-0.12.1.yaml"
    benchmark:
        os.path.join(config["OUT_ROOT"], "benchmarks", "run_fastqc.{sample}.txt")
    threads: 12
    # Unlike the demux, a repeat run of this rule can succeed: it is idempotent
    # and its usual failures are transient (node or filesystem), so the
    # attempt-scaled resources below do real work.
    retries: 3
    resources:
        mem_gb=lambda wildcards, attempt: 16 + ((attempt - 1) * 32),
        time=lambda wildcards, attempt: 180 + ((attempt - 1) * 720),
        partition="compute,hugemem,vgpu",
    shell:
        """ 
        
        mkdir -p {fastqc_out_root}

        fastqc -t {threads} -o {fastqc_out_root} {input.fastq} > {log} 2>&1

        success_landmark={output.zip}

        if [ ! -f $success_landmark ]
        then
            echo "error: fastqc  did not generate the expected output file {output.zip}. "
            exit 1
        else
            exit 0
        fi

        """


def get_fastq_reports(wildcards, extension=".fastq.gz"):
    directory = checkpoints.run_bclconvert.get().output["bclconvert_out"]
    files = [
        os.path.basename(f)
        for f in os.listdir(directory)
        if f.endswith(extension) and "Undetermined" not in f
    ]
    basenames = [f.replace(extension, "") for f in files]
    return expand(
        os.path.join(config["OUT_ROOT"], "SampleSheet", "fastqc_run", "fastqc", "{sample}_fastqc.zip"), 
        sample=basenames
    )


multiqc_data_dir_path = os.path.join(
    config["OUT_ROOT"], "SampleSheet", "multiqc", "multiqc_data"
)

localrules: run_multiqc
rule run_multiqc:
    input:
        fastqc_reports=lambda wildcards: get_fastq_reports(wildcards),
        bclconvert_reports=os.path.join(config["OUT_ROOT"], "SampleSheet", "bclconvert", "Reports"),
    output:
        report=os.path.join(
            config["OUT_ROOT"],
            "SampleSheet",
            "multiqc",
            f"{config['RUN']}.multiqc.html",
        ),
    log:
        os.path.join(config["OUT_ROOT"], "logs", "run_multiqc.log"),
    conda:
        "envs/multiqc-1.17.yaml"
    benchmark:
        os.path.join(config["OUT_ROOT"], "benchmarks", "run_multiqc.txt")
    # A localrule on the submit host, and the last step after hours of work.
    # One retry, because a MultiQC failure is usually deterministic and a
    # second identical attempt rarely helps.
    retries: 1
    threads: 2
    params:
        multiqc_config=config["multiqc_config"],
    shell:
        """
        multiqc --interactive --outdir {multiqc_data_dir_path} --filename {output.report} --force -c {params.multiqc_config} --data-dir --data-format tsv {input.bclconvert_reports} {input.fastqc_reports} > {log} 2>&1
        """
