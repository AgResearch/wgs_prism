# 2026 Benjamin J Perry
# MIT License
# Copyright (c) 2026 Benjamin J Perry
# Version: 1.0
# Maintainer: Benjamin J Perry
# Email: ben.perry@agresearch.co.nz


configfile: "config/pipeline_config.yaml"


import os

# config dictionary values to be defined on running snakemake with --config flag

RUN = config["RUN"]
LANE = config["LANE"]

# MGI emits one multiplexed fastq per lane, with the barcodes still on the reads.
fastq_in = os.path.join(config["IN_ROOT"], RUN, LANE, f"{RUN}_{LANE}_read.fq.gz")
run_info_in = os.path.join(
    config["IN_ROOT"], RUN, "Metrics", f"{RUN}_RunInfoMetrics.csv"
)

demux_dir = os.path.join(config["OUT_ROOT"], "splitbarcode")
fastqc_dir = os.path.join(config["OUT_ROOT"], "fastqc")

# splitBarcode names its outputs <RUN>_<LANE>_<Sample_ID>.fq.gz from the read
# header, so FastQC reports (and therefore MultiQC samples) carry this prefix.
sample_prefix = f"{RUN}_{LANE}_"

# splitBarcode is not packaged for bioconda; it comes from the module system.
SPLITBARCODE_ROOT = "/agr/persist/apps/eri_rocky8/software/MGI-splitBarcode/2.0.0-4"


rule targets:
    input:
        os.path.join(config["OUT_ROOT"], "multiqc", f"{RUN}.multiqc.html"),


localrules:
    prepare_barcodes,
    barcodestat_to_mqc,
    run_multiqc,


rule prepare_barcodes:
    input:
        sample_sheet=config["SAMPLE_SHEET"],
        run_info=run_info_in,
        fastq=fastq_in,
    output:
        barcodes=os.path.join(config["OUT_ROOT"], f"{RUN}.barcode"),
        params=os.path.join(config["OUT_ROOT"], "splitbarcode_params.json"),
    log:
        os.path.join(config["OUT_ROOT"], "logs", "prepare_barcodes.log"),
    benchmark:
        os.path.join(config["OUT_ROOT"], "benchmarks", "prepare_barcodes.txt")
    params:
        mismatch=config["barcode_mismatch"],
        lane=LANE,
    shell:
        """
        python3 workflow/scripts/prepare_barcodes.py \
            --sample-sheet {input.sample_sheet} \
            --lane {params.lane} \
            --run-info {input.run_info} \
            --fastq {input.fastq} \
            --out-barcodes {output.barcodes} \
            --out-params {output.params} \
            --mismatch {params.mismatch} > {log} 2>&1
        """


# splitBarcode writes its own log to ./log/ relative to the working directory,
# not to -o. Running it from here therefore keeps that log outside the
# directory() output (so it survives Snakemake deleting that output on failure)
# and keeps it out of the checkout, which is the working directory otherwise.
# Note demux_dir itself is only created as a side effect of preparing the
# nested BarcodeStat.txt output below (Snakemake pre-creates a directory()
# output's parent, not the directory itself) — splitBarcode writes happily into
# the resulting empty directory.
splitbarcode_logs_persist = os.path.join(config["OUT_ROOT"], "logs", "splitbarcode_Logs")


rule run_splitbarcode:
    input:
        fastq=fastq_in,
        barcodes=os.path.join(config["OUT_ROOT"], f"{RUN}.barcode"),
        params=os.path.join(config["OUT_ROOT"], "splitbarcode_params.json"),
    output:
        demux=directory(demux_dir),
        # Declared explicitly so barcodestat_to_mqc can depend on it directly.
        barcodestat=os.path.join(demux_dir, "BarcodeStat.txt"),
    log:
        os.path.join(config["OUT_ROOT"], "logs", "run_splitbarcode.log"),
    benchmark:
        os.path.join(config["OUT_ROOT"], "benchmarks", "run_splitbarcode.txt")
    threads: 16
    resources:
        mem_gb=lambda wildcards, attempt: 120 + ((attempt - 1) * 40),
        time=lambda wildcards, attempt: 480 + ((attempt - 1) * 480),
        partition="compute,hugemem,vgpu",
    shell:
        """
        export PATH={SPLITBARCODE_ROOT}/bin:$PATH
        export LD_LIBRARY_PATH={SPLITBARCODE_ROOT}/lib:$LD_LIBRARY_PATH

        mkdir -p {output.demux}
        mkdir -p {splitbarcode_logs_persist}

        # splitBarcode is run from the log dir below, so give it absolute paths.
        sb_barcodes=$(realpath {input.barcodes})
        sb_fastq=$(realpath {input.fastq})
        sb_demux=$(realpath {output.demux})

        # Build -b/-r from the validated params so the geometry lives in one
        # place. Read at runtime rather than via a params: lambda, which
        # Snakemake may evaluate before prepare_barcodes has written the file.
        sb_args=$(python3 -c "
import json
p = json.load(open('{input.params}'))
print(' '.join('-b %d %d %d' % tuple(b) for b in p['b']) + (' -r' if p['reverse'] else ''))
")

        echo "splitBarcode arguments: $sb_args" > {log}
        echo >> {log}

        # Capture exit status ourselves (Snakemake's `set -e` would otherwise
        # abort before the log below is appended on failure).
        sb_status=0
        ( cd {splitbarcode_logs_persist} && splitBarcode -B $sb_barcodes \
            -1 $sb_fastq \
            $sb_args \
            -t {threads} -m 100 \
            -o $sb_demux ) >> {log} 2>&1 || sb_status=$?

        # splitBarcode's own log/ is already outside the directory() output, so
        # it survives failure; surface it in the rule log too. The `ls` guard
        # avoids tripping `set -e` if it never appeared.
        if ls {splitbarcode_logs_persist}/log/* >/dev/null 2>&1; then
            echo >> {log}
            echo "===== splitBarcode log/ (preserved in {splitbarcode_logs_persist}/log) =====" >> {log}
            cat {splitbarcode_logs_persist}/log/* >> {log} 2>/dev/null || true
        fi

        exit $sb_status
        """


rule run_fastqc:
    input:
        demux=demux_dir,
    output:
        fastqc=directory(fastqc_dir),
    log:
        os.path.join(config["OUT_ROOT"], "logs", "run_fastqc.log"),
    conda:
        "envs/fastqc-0.12.1.yaml"
    benchmark:
        os.path.join(config["OUT_ROOT"], "benchmarks", "run_fastqc.txt")
    threads: 32
    resources:
        mem_gb=lambda wildcards, attempt: 32 + ((attempt - 1) * 32),
        time=lambda wildcards, attempt: 120 + ((attempt - 1) * 240),
        partition="compute,hugemem,vgpu",
    shell:
        """
        mkdir -p {output.fastqc}

        # FastQC does not recurse a directory, so glob. Exclude undecoded (the
        # MGI equivalent of Illumina's Undetermined), and pipe via xargs rather
        # than expanding several thousand paths onto one command line.
        find {input.demux} -name '*.fq.gz' ! -name '*undecoded*' \
            | xargs fastqc -t {threads} --noextract -o {output.fastqc} > {log} 2>&1

        # Every input must have produced a report. Counting rather than probing
        # for one landmark file also catches a partial run. Do not pipe a long
        # `ls` into `head` here: it takes SIGPIPE, and `set -o pipefail` then
        # fails the rule after a perfectly good FastQC run.
        n_expected=$(find {input.demux} -name '*.fq.gz' ! -name '*undecoded*' | wc -l)
        n_reports=$(find {output.fastqc} -name '*_fastqc.zip' | wc -l)

        if [ "$n_reports" -ne "$n_expected" ]
        then
            echo "error: fastqc produced $n_reports reports in {output.fastqc}, expected $n_expected."
            exit 1
        else
            exit 0
        fi
        """


rule barcodestat_to_mqc:
    input:
        barcodestat=os.path.join(demux_dir, "BarcodeStat.txt"),
    output:
        mqc=os.path.join(config["OUT_ROOT"], "multiqc", "splitbarcode_mqc.tsv"),
    log:
        os.path.join(config["OUT_ROOT"], "logs", "barcodestat_to_mqc.log"),
    params:
        sample_prefix=sample_prefix,
    shell:
        """
        python3 workflow/scripts/barcodestat_to_mqc.py \
            --barcodestat {input.barcodestat} \
            --out {output.mqc} \
            --sample-prefix {params.sample_prefix} > {log} 2>&1
        """


multiqc_data_dir_path = os.path.join(config["OUT_ROOT"], "multiqc", "multiqc_data")


rule run_multiqc:
    input:
        fastqc=fastqc_dir,
        splitbarcode_mqc=os.path.join(
            config["OUT_ROOT"], "multiqc", "splitbarcode_mqc.tsv"
        ),
    output:
        report=os.path.join(config["OUT_ROOT"], "multiqc", f"{RUN}.multiqc.html"),
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
        multiqc --interactive --outdir {multiqc_data_dir_path} --filename {output.report} --force -c {params.multiqc_config} --data-dir --data-format tsv {input.splitbarcode_mqc} {input.fastqc} > {log} 2>&1
        """
