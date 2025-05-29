# 2024 Benjamin J Perry
# MIT License
# Copyright (c) 2024 Benjamin J Perry
# Version: 1.0
# Maintainer: Benjamin J Perry
# Email: ben.perry@agresearch.co.nz


configfile: "config/pipeline_config.yaml"


import os
import pandas as pd


# Global variables
run_name = config["RUN"]
run_in_path = os.path.join(config["IN_ROOT"], config["RUN"])
bclconvert_out_root = os.path.join(config["OUT_ROOT"])
bclconvert_out_path = os.path.join(bclconvert_out_root, "SampleSheet/bclconvert")


# config dictionary values to be defined on running snakemake with --config flag
# Note: until fully migrated, singularity mounting doesn't work for this as there are too many symlinks. The actual run_in_path will need to be updated to /mnt/gpfs/persist/projects/2024_illumina_sequencing_d/active
sample_sheet_path = os.path.join(bclconvert_out_root, "SampleSheet.csv")
bclconvert_out_path = os.path.join(bclconvert_out_root, "SampleSheet/bclconvert")
top_unknown_path = os.path.join(
    bclconvert_out_root, "SampleSheet/bclconvert/Reports/Top_Unknown_Barcodes.csv"
)
fastq_complete_path = os.path.join(
    bclconvert_out_root, "SampleSheet/bclconvert/Logs/FastqComplete.txt"
)
bclconvert_log = os.path.join(bclconvert_out_root, "logs/1_run_bclconvert.log")
bclconvert_benchmark = os.path.join(
    bclconvert_out_root, "benchmarks/run_bclconvert.txt"
)


# config dictionary values to be defined on running snakemake with --config flag
fastqc_in_root = os.path.join(config["OUT_ROOT"], "SampleSheet/bclconvert")
fastqc_in_samples = os.path.join(
    config["OUT_ROOT"], "SampleSheet/bclconvert/{sample}.fastq.gz"
)
fastqc_out_root = os.path.join(config["OUT_ROOT"], "SampleSheet/fastqc_run/fastqc")
fastqc_out_samples_zips = os.path.join(
    config["OUT_ROOT"], "SampleSheet/fastqc_run/fastqc/{sample}_fastqc.zip"
)
fastqc_out_samples_htmls = os.path.join(
    config["OUT_ROOT"], "SampleSheet/fastqc_run/fastqc/{sample}_fastqc.html"
)
fastqc_log = os.path.join(config["OUT_ROOT"], "logs/2.1_run_fastqc.{sample}.log")
fastqc_benchmark = os.path.join(
    config["OUT_ROOT"], "benchmarks/run_fastqc.{sample}.txt"
)


# logs and reports for multiQC
bclconvert_reports_dir = os.path.join(
    config["OUT_ROOT"], "SampleSheet/bclconvert/Reports"
)
fastqc_reports_dir = os.path.join(config["OUT_ROOT"], "SampleSheet/fastqc_run/fastqc")
multiqc_report_file = run_name + ".multiqc.html"
multiqc_report_path = os.path.join(
    config["OUT_ROOT"], "SampleSheet", "multiqc", multiqc_report_file
)
multiqc_data_dir = "multiqc_data"
multiqc_data_dir_path = os.path.join(
    config["OUT_ROOT"], "SampleSheet", "multiqc", multiqc_data_dir
)
multiqc_log = "logs/3.0.0_run_multiqc.log"
multiqc_log_path = os.path.join(config["OUT_ROOT"], multiqc_log)
multiqc_benchmark = "benchmarks/run_multiqc.txt"
multiqc_benchmark_path = os.path.join(config["OUT_ROOT"], multiqc_benchmark)


rule targets:
    input:
        multiqc_report_path,


checkpoint run_bclconvert:
    input:
        run_in=run_in_path,
        sample_sheet=sample_sheet_path,
    output:
        bclconvert_out=directory(bclconvert_out_path),
        fastq_complete=fastq_complete_path,
    benchmark:
        bclconvert_benchmark
    threads: 36
    resources:
        mem_gb=lambda wildcards, attempt: 128 + ((attempt - 1) * 32),
        time=lambda wildcards, attempt: 480 + ((attempt - 1) * 120),
    shell:
        """
        export PATH=/agr/persist/apps/src/b/BCL-Convert:$PATH # For ERI

        # report version 
        echo "bcl-convert version in use:"
        bcl-convert -V 

        echo

        bcl-convert --bcl-input-directory {input.run_in} --sample-sheet {input.sample_sheet} --output-directory {output.bclconvert_out}

        """


rule run_fastqc:
    input:
        fastq=fastqc_in_samples,
    output:
        zip=fastqc_out_samples_zips,
        html=fastqc_out_samples_htmls,
    log:
        fastqc_log,
    conda:
        "envs/fastqc-0.12.1.yaml"
    benchmark:
        fastqc_benchmark
    threads: 12
    resources:
        mem_gb=lambda wildcards, attempt: 8 + ((attempt - 1) * 32),
        time=lambda wildcards, attempt: 120 + ((attempt - 1) * 120),
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


def get_fastq(wildcards, extension=".fastq.gz"):
    directory = checkpoints.run_bclconvert.get().output[0]
    files = [
        subpath(f, basename=True)
        for f in os.listdir(directory)
        if f.endswith(extension)
    ]
    basenames = [f.replace(extension, "") for f in files]
    return expand(
        os.path.join("results", library, "01_cutadapt/{samples}.fastq.gz"),
        samples=basenames,
    )


rule run_multiqc:
    input:
        expand(fastqc_out_samples_zips, sample=get_fastq),  # input function for fastq results
        bclconvert_in=bclconvert_reports_dir,
    output:
        report=multiqc_report_path,
    log:
        multiqc_log_path,
    conda:
        "envs/multiqc-1.17.yaml"
    benchmark:
        multiqc_benchmark_path
    threads: 2
    resources:
        mem_gb=lambda wildcards, attempt: 16 + ((attempt - 1) * 8),
        time=lambda wildcards, attempt: 30 + ((attempt - 1) * 30),
    params:
        multiqc_config=config["multiqc_config"],
    shell:
        """

        multiqc --interactive --outdir {multiqc_data_dir_path} --filename {output.report} --force -c {params.multiqc_config} --data-dir --data-format tsv {input.bclconvert_in} {input.fastqc_in} > {log} 2>&1
        """
