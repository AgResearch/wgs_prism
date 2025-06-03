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

rule targets:
    input:
        os.path.join(
            config["OUT_ROOT"],
            "SampleSheet",
            "multiqc",
            f"{config['RUN']}.multiqc.html",
        ),


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
    resources:
        mem_gb=lambda wildcards, attempt: 96 + ((attempt - 1) * 24),
        time=lambda wildcards, attempt: 480 + ((attempt - 1) * 120),
        partition="compute,hugemem,vgpu",
    shell:
        """
        export PATH=/agr/persist/apps/src/b/BCL-Convert:$PATH

        # report version 
        echo "bcl-convert version in use:"
        bcl-convert -V 

        echo

        bcl-convert --force --bcl-input-directory {input.run_in} --sample-sheet {input.sample_sheet} --output-directory {output.bclconvert_out} > {log} 2>&1

        # Need to scrape the logs into the log file for the failure case as snakemake cleans up the directory with the logs
        cat {output.bclconvert_out}/Logs/*.log > {log}
        
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
    resources:
        mem_gb=lambda wildcards, attempt: 8 + ((attempt - 1) * 32),
        time=lambda wildcards, attempt: 120 + ((attempt - 1) * 120),
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
        if f.endswith(extension)
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
    threads: 2
    params:
        multiqc_config=config["multiqc_config"],
    shell:
        """
        multiqc --interactive --outdir {multiqc_data_dir_path} --filename {output.report} --force -c {params.multiqc_config} --data-dir --data-format tsv {input.bclconvert_reports} {input.fastqc_reports} > {log} 2>&1
        """
