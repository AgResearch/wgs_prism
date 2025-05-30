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


rule targets:
    input:
        os.path.join(
            config["OUT_ROOT"],
            "SampleSheet",
            "multiqc",
            (config["RUN"] + ".multiqc.html"),
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
    benchmark:
        os.path.join(config["OUT_ROOT"], "benchmarks", "run_bclconvert.txt")
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
    return basenames


multiqc_data_dir_path = os.path.join(
    config["OUT_ROOT"], "SampleSheet", "multiqc", "multiqc_data"
)


rule run_multiqc:
    input:
        fastqc_reports=lambda wildcards: expand(
            rules.run_fastqc.output.zip, sample=get_fastq(wildcards)
        ),
        bclconvert_reports=checkpoints.run_bclconvert.get().output["reports"],
    output:
        report=os.path.join(
            config["OUT_ROOT"],
            "SampleSheet",
            "multiqc",
            (config["RUN"] + ".multiqc.html"),
        ),
    log:
        os.path.join(config["OUT_ROOT"], "logs", "run_multiqc.log"),
    conda:
        "envs/multiqc-1.17.yaml"
    benchmark:
        os.path.join(config["OUT_ROOT"], "benchmarks", "run_multiqc.txt")
    threads: 2
    resources:
        mem_gb=lambda wildcards, attempt: 16 + ((attempt - 1) * 8),
        time=lambda wildcards, attempt: 30 + ((attempt - 1) * 30),
    params:
        multiqc_config=config["multiqc_config"],
    shell:
        """
        multiqc --interactive --outdir {multiqc_data_dir_path} --filename {output.report} --force -c {params.multiqc_config} --data-dir --data-format tsv {input.bclconvert_reports} {input.fastqc_reports} > {log} 2>&1
        """
