# 2024 Benjamin J Perry
# MIT License
# Copyright (c) 2024 Benjamin J Perry
# Version: 1.0
# Maintainer: Benjamin J Perry
# Email: ben.perry@agresearch.co.nz

configfile: 'config/pipeline_config.yaml'

import os
import pandas as pd


# onstart:
#     print(f"Working directory: {os.getcwd()}")
#     print("TOOLS: ")
#     os.system('echo "  bash: $(which bash)"')
#     os.system('echo "  PYTHON: $(which python)"')
#     os.system('echo "  CONDA: $(which conda)"')
#     os.system('echo "  SNAKEMAKE: $(which snakemake)"')
#     print(f"Env TMPDIR = {os.environ.get('TMPDIR', '<n/a>')}")
#     os.system('echo "  PYTHON VERSION: $(python --version)"')
#     os.system('echo "  CONDA VERSION: $(conda --version)"')


# wildcard_constraints:
#     samples="\w+"

# Global variables
# config dictionary values to be defined on running snakemake with --config flag
# Note: until fully migrated, singularity mounting doesn't work for this as there are too many symlinks. The actual run_in_path will need to be updated to /mnt/gpfs/persist/projects/2024_illumina_sequencing_d/active
run_in_path = os.path.join(config["IN_ROOT"], config["RUN"])
bclconvert_out_root = os.path.join(config["OUT_ROOT"])

sample_sheet_path = os.path.join(bclconvert_out_root, "SampleSheet.csv")
bclconvert_out_path = os.path.join(bclconvert_out_root, "SampleSheet/bclconvert")
top_unknown_path = os.path.join(bclconvert_out_root, "SampleSheet/bclconvert/Reports/Top_Unknown_Barcodes.csv")
fastq_complete_path = os.path.join(bclconvert_out_root, "SampleSheet/bclconvert/Logs/FastqComplete.txt")

bclconvert_log = os.path.join(bclconvert_out_root, "logs/1_run_bclconvert.log")
bclconvert_benchmark = os.path.join(bclconvert_out_root, "benchmarks/run_bclconvert.txt")


rule targets:
    input:
        fastq_complete_path,
        # top_unknown_path


rule run_bclconvert:
    input:
        run_in = run_in_path,
        sample_sheet = sample_sheet_path,
    output:
        bclconvert_out = directory(bclconvert_out_path),
        fastq_complete = fastq_complete_path,
        # top_unknown = top_unknown_path
    log:
        bclconvert_log
    # singularity:
    #     "docker://nfcore/bclconvert:3.9.3"
    benchmark:
        bclconvert_benchmark
    threads: 36
    resources:
        mem_gb = lambda wildcards, attempt: 128 + ((attempt - 1) * 32),
        time = lambda wildcards, attempt: 480 + ((attempt - 1) * 120),
    shell:
        """
        export PATH=/agr/persist/apps/src/b/BCL-Convert:$PATH # For ERI
        # run bcl-convert
        # report version 
        touch {bclconvert_log}
        echo "bcl-convert version in use:"

        bcl-convert -V 

        echo

        bcl-convert --force --bcl-input-directory {input.run_in} --sample-sheet {input.sample_sheet} --output-directory {output.bclconvert_out} 2>&1 | tee -a {log} 

        cat {bclconvert_out_path}/Logs/*log >> {log}

        """
        
