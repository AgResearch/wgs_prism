# 2024 Benjamin J Perry
# MIT License
# Copyright (c) 2024 Benjamin J Perry
# Version: 1.0
# Maintainer: Benjamin J Perry
# Email: ben.perry@agresearch.co.nz
# configfile: "config/config.yaml"


import os
configfile: 'config/pipeline_config.yaml'

# eg. /dataset/2024_illumina_sequencing_d/active/240621_A01439_0277_BHVHYHDMXY
RUN = config["RUN"]
run_in_path = os.path.join(config["IN_ROOT"], config["RUN"]) 

# eg. /dataset/2024_illumina_sequencing_d/scratch/postprocessing/illumina/novaseq/240621_A01439_0277_BHVHYHDMXY/SampleSheet.csv
sample_sheet_path = os.path.join(config["OUT_ROOT"], "SampleSheet.csv") 

# eg. /dataset/2024_illumina_sequencing_d/scratch/postprocessing/illumina/novaseq/240621_A01439_0277_BHVHYHDMXY/SampleSheet
OUT_ROOT = os.path.join(config["OUT_ROOT"], "SampleSheet") 

# eg. /dataset/2024_illumina_sequencing_d/scratch/postprocessing/illumina/novaseq/240621_A01439_0277_BHVHYHDMXY/logs
LOG_ROOT = os.path.join(config["OUT_ROOT"], "logs")

# eg. /dataset/2024_illumina_sequencing_d/scratch/postprocessing/illumina/novaseq/240621_A01439_0277_BHVHYHDMXY/benchmarks
BENCHMARK_ROOT = os.path.join(config["OUT_ROOT"], "benchmarks")

bclconvert_out_path = os.path.join(OUT_ROOT, "bclconvert")
bclconvert_reports = os.path.join(OUT_ROOT, "bclconvert/Reports")


multiqc_report_file = RUN + ".multiqc.html"
multiqc_report_path = os.path.join(OUT_ROOT, multiqc_report_file)


onstart:
    print(f"Working directory: {os.getcwd()}")
    print("TOOLS: ")
    os.system('echo "  bash: $(which bash)"')
    os.system('echo "  PYTHON: $(which python)"')
    os.system('echo "  CONDA: $(which conda)"')
    os.system('echo "  SNAKEMAKE: $(which snakemake)"')
    print(f"Env TMPDIR = {os.environ.get('TMPDIR', '<n/a>')}")
    os.system('echo "  PYTHON VERSION: $(python --version)"')
    os.system('echo "  CONDA VERSION: $(conda --version)"')


wildcard_constraints: sample = "(?!Undetermined).+"

# FIDs is used in expands for multiQC report
fastqc_in_root = os.path.join(OUT_ROOT, "bclconvert")
(FIDs,) = glob_wildcards( os.path.join(fastqc_in_root, "{sample,(?!Undetermined).*}_R1_001.fastq.gz") )


rule default:
    input:
        multiqc_report_path,


rule downsample_for_QC:
    input:
        read1 = os.path.join(OUT_ROOT, "bclconvert", "{samples}_R1_001.fastq.gz"),
        read2 = os.path.join(OUT_ROOT, "bclconvert", "{samples}_R2_001.fastq.gz"),
    output:
        downed1 = os.path.join(OUT_ROOT, "00_downsample", "{samples}.DS.R1.fastq.gz"),
        downed2 = os.path.join(OUT_ROOT, "00_downsample", "{samples}.DS.R2.fastq.gz"),
    log:
        os.path.join(LOG_ROOT, "seqkit.{samples}.log")
    benchmark:
        os.path.join(BENCHMARK_ROOT, "seqkit.{samples}.txt")
    conda:
        "seqtk"
    threads: 2
    resources:
        mem_gb = lambda wildcards, attempt: 2 + ((attempt - 1) * 4),
        time = lambda wildcards, attempt: 6 + ((attempt - 1) * 6),
    shell:
        "seqtk sample -s1953 {input.read1} 1000000 | gzip > {output.downed1} 2> {log} "
        "&& "
        "seqtk sample -s1953 {input.read2} 1000000 | gzip > {output.downed2} 2> {log} "


rule bbduk_read_trim:
    input:
        downed1 = os.path.join(OUT_ROOT, "00_downsample", "{samples}.DS.R1.fastq.gz"),
        downed2 = os.path.join(OUT_ROOT, "00_downsample", "{samples}.DS.R2.fastq.gz"),
        adapters = "resources/adapters.fa",
    output:
        bbdukRead1 = os.path.join(OUT_ROOT, "01_readMasking", "{samples}_R1_bbduk.fastq.gz"),
        bbdukRead2 = os.path.join(OUT_ROOT, "01_readMasking", "{samples}_R2_bbduk.fastq.gz"),
    log:
        os.path.join(LOG_ROOT, "bbduk.{samples}.log")
    benchmark:
        os.path.join(BENCHMARK_ROOT, "bbduk.{samples}.txt")
    conda:
        "bbduk"
    threads: 8
    resources: #TODO Update
        mem_gb = lambda wildcards, attempt: 6 + ((attempt - 1) * 6),
        time = lambda wildcards, attempt: 20 + ((attempt - 1) * 15),
    shell:
        "bbduk.sh "
        "threads={threads} "
        "in1={input.downed1} "
        "in2={input.downed2} "
        "ref={input.adapters} "
        "ktrim=rl "
        "k=19 mink=9 hdist=1 "
        "rcomp=t " 
        "trimpolyg=2 "
        "trimpolya=10 "
        "qtrim=f "
        "minlen=10 "
        "out1={output.bbdukRead1} "
        "out2={output.bbdukRead2} "
        "2>&1 | tee {log} "

rule bowtie2_SILVA_alignment_read1:
    input:
        bbdukRead1 = os.path.join(OUT_ROOT, "01_readMasking", "{samples}_R1_bbduk.fastq.gz"),
    output:
        bowtie2_R1 = os.path.join(OUT_ROOT, "02_SILVA", "{samples}.DS.R1.bowtie2.log"),
#        silva_R1 = "results/02_SILVA/{samples}_R1_bbduk_silva.fastq",
    benchmark:
        os.path.join(BENCHMARK_ROOT, "bowtie2_SILVA_alignment_read1.{samples}.txt")
    conda:
        "bowtie2-2.5.1"
    threads: 8
    resources:
        mem_gb = lambda wildcards, attempt: 24 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 10 + ((attempt - 1) * 20),
    shell:
        "bowtie2 "
        "-p {threads} "
        "-x /datasets/2024-silva-rrna/SILVA138.1 " # TODO parameterize in config
#        "--un {output.silva_R1} "
        "-U {input.bbdukRead1} "
        "1> /dev/null "
        "2> {output.bowtie2_R1} "


rule bowtie2_SILVA_alignment_read2:
    input:
        bbdukRead2 = os.path.join(OUT_ROOT, "01_readMasking", "{samples}_R2_bbduk.fastq.gz"),
    output:
        bowtie2_R2 = os.path.join(OUT_ROOT, "02_SILVA", "{samples}.DS.R2.bowtie2.log"),
    benchmark:
        os.path.join(BENCHMARK_ROOT, "bowtie2_SILVA_alignment_read2.{samples}.txt")
    conda:
        "bowtie2-2.5.1"
    threads: 8
    resources:
        mem_gb = lambda wildcards, attempt: 24 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 10 + ((attempt - 1) * 20),
    shell:
        "bowtie2 "
        "-p {threads} "
        "-x /datasets/2024-silva-rrna/SILVA138.1 "
#        "--un {output.silva_R2} "
        "-U {input.bbdukRead2} "
        "1> /dev/null "
        "2> {output.bowtie2_R2} "

kraken2_index = config['kraken2_index']

rule kraken2_read_composition_read1:
    input:
        bbdukRead1 = os.path.join(OUT_ROOT, "01_readMasking", "{samples}_R1_bbduk.fastq.gz"),
    output:
        k2Output = temp(os.path.join(OUT_ROOT, "03_kraken2", "{samples}.DS.R1.nt.k2")),
        k2Report_R1 = os.path.join(OUT_ROOT, "03_kraken2", "{samples}.DS.R1.nt.report.kraken2"),
    log:
        os.path.join(LOG_ROOT, "kraken2_read1_composition.{samples}.log"),
    benchmark:
        os.path.join(BENCHMARK_ROOT, "kraken2_read1_composition.{samples}.txt"),
    conda:
        "kraken2"
    threads: 16
    resources:
        mem_gb = lambda wildcards, attempt: 714 + ((attempt - 1) * 20),
        time = lambda wildcards, attempt: 35 + ((attempt - 1) * 40),
        partition = "hugemem",
    shell:
        "kraken2 "
        "--use-names "
        "--db {kraken2_index} " 
        "-t {threads} "
        "--report {output.k2Report_R1} "
        "--report-minimizer-data "
        "--output {output.k2Output} "
        "{input.bbdukRead1} "
        "2>&1 | tee {log} "


rule kraken2_read_composition_read2:
    input:
        bbdukRead2 = os.path.join(OUT_ROOT, "01_readMasking", "{samples}_R2_bbduk.fastq.gz"),
    output:
        k2Output = temp(os.path.join(OUT_ROOT, "03_kraken2", "{samples}.DS.R2.nt.k2")),
        k2Report_R2 = os.path.join(OUT_ROOT, "03_kraken2", "{samples}.DS.R2.nt.report.kraken2"),
    log:
        os.path.join(LOG_ROOT, "kraken2_read2_composition.{samples}.log"),
    benchmark:
        os.path.join(BENCHMARK_ROOT, "kraken2_read2_composition.{samples}.txt"),
    conda:
        "kraken2"
    threads: 16
    resources:
        mem_gb = lambda wildcards, attempt: 714 + ((attempt - 1) * 20),
        time = lambda wildcards, attempt: 35 + ((attempt - 1) * 40),
        partition = "hugemem",
    shell:
        "kraken2 "
        "--use-names "
        "--db {kraken2_index} " 
        "-t {threads} "
        "--report {output.k2Report_R2} "
        "--report-minimizer-data "
        "--output {output.k2Output} "
        " {input.bbdukRead2} "
        "2>&1 | tee {log} "


fastqc_out_root = os.path.join(OUT_ROOT, "fastqc_run", "fastqc")

rule fastqc_read1:
    input:
        read1 = os.path.join(OUT_ROOT, "bclconvert", "{samples}_R1_001.fastq.gz"),
    output:
        html = os.path.join(OUT_ROOT, "fastqc_run", "fastqc", "{samples}_R1_001_fastqc.html"),
        zip = os.path.join(OUT_ROOT, "fastqc_run", "fastqc", "{samples}_R1_001_fastqc.zip"),
    conda:
        "fastqc-0.12.1"
    threads: 4
    resources:
        mem_gb = lambda wildcards, attempt: 6 + ((attempt - 1) * 2),
        time = lambda wildcards, attempt: 20 + ((attempt - 1) * 15),
    shell:
        "fastqc "
        "-o {fastqc_out_root} "
        "-q "
        "-t {threads} "
        "{input.read1}"


rule fastqc_read2:
    input:
        read2 = os.path.join(OUT_ROOT, "bclconvert", "{samples}_R2_001.fastq.gz"),
    output:
        html = os.path.join(OUT_ROOT, "fastqc_run", "fastqc", "{samples}_R2_001_fastqc.html"),
        zip = os.path.join(OUT_ROOT, "fastqc_run", "fastqc", "{samples}_R2_001_fastqc.zip"),
    conda:
        "fastqc-0.12.1"
    threads: 4
    resources:
        mem_gb = lambda wildcards, attempt: 6 + ((attempt - 1) * 2),
        time = lambda wildcards, attempt: 20 + ((attempt - 1) * 15),
    shell:
        "fastqc "
        "-o {fastqc_out_root} "
        "-q "
        "-t {threads} "
        "{input.read2}"


rule fastqc_filtered_read1: 
    input:
        bbdukRead1 = os.path.join(OUT_ROOT, "01_readMasking", "{samples}_R1_bbduk.fastq.gz"),
    output:
        html = os.path.join(OUT_ROOT, "fastqc_run", "fastqc", "{samples}_R1_bbduk_fastqc.html"),
        zip = os.path.join(OUT_ROOT, "fastqc_run", "fastqc", "{samples}_R1_bbduk_fastqc.zip"),
    conda:
        "fastqc-0.12.1"
    threads: 4
    resources:
        mem_gb = lambda wildcards, attempt: 6 + ((attempt - 1) * 2),
        time = lambda wildcards, attempt: 20 + ((attempt - 1) * 15),
    shell:
        "fastqc "
        "-o {fastqc_out_root} "
        "-q "
        "-t {threads} "
        "{input.bbdukRead1}"


rule fastqc_filtered_read2: 
    input:
        bbdukRead2 = os.path.join(OUT_ROOT, "01_readMasking", "{samples}_R2_bbduk.fastq.gz"),
    output:
        html = os.path.join(OUT_ROOT, "fastqc_run", "fastqc", "{samples}_R2_bbduk_fastqc.html"),
        zip = os.path.join(OUT_ROOT, "fastqc_run", "fastqc", "{samples}_R2_bbduk_fastqc.zip"),
    conda:
        "fastqc-0.12.1"
    threads: 4
    resources:
        mem_gb = lambda wildcards, attempt: 6 + ((attempt - 1) * 2),
        time = lambda wildcards, attempt: 20 + ((attempt - 1) * 15),
    shell:
        "fastqc "
        "-o {fastqc_out_root} "
        "-q "
        "-t {threads} "
        "{input.bbdukRead2}"

human_genome_index = config['human_genome_index']

rule genome_alignment_check_R1:
    input:
        bbdukRead1 = os.path.join(OUT_ROOT, "01_readMasking", "{samples}_R1_bbduk.fastq.gz"),
    output:
        bowtie2_genome = os.path.join(OUT_ROOT, "02_REF", "{samples}.DS.genome_alignment.bowtie2.R1.log"),
    benchmark:
        os.path.join(BENCHMARK_ROOT, "genome_alignment_check.R1.{samples}.txt"),
    conda:
        "bowtie2-2.5.1"
    threads: 8
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 10 + ((attempt - 1) * 30),
    shell:
        "bowtie2 "
        "-p {threads} "
        "-x {human_genome_index} "
        "-U {input.bbdukRead1} "
        "1> /dev/null "
        "2> {output.bowtie2_genome} "


rule genome_alignment_check_R2:
    input:
        bbdukRead2 = os.path.join(OUT_ROOT, "01_readMasking", "{samples}_R2_bbduk.fastq.gz"),
    output:
        bowtie2_genome = os.path.join(OUT_ROOT, "02_REF", "{samples}.DS.genome_alignment.bowtie2.R2.log"),
    benchmark:
        os.path.join(BENCHMARK_ROOT, "genome_alignment_check.R2.{samples}.txt"),
    conda:
        "bowtie2-2.5.1"
    threads: 8
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 10 + ((attempt - 1) * 30),
    shell:
        "bowtie2 "
        "-p {threads} "
        "-x {human_genome_index} "
        "-U {input.bbdukRead2} "
        "1> /dev/null "
        "2> {output.bowtie2_genome} "


rule multiQC_report:
    input:
        bclconvert_runinfo = os.path.join(bclconvert_reports, "RunInfo.xml") ,
        bclconvert_demux_stats = os.path.join(bclconvert_reports, "Demultiplex_Stats.csv"),
        bclconvert_qual_metrics = os.path.join(bclconvert_reports, "Quality_Metrics.csv"),
        bclconvert_adpt_metrics = os.path.join(bclconvert_reports, "Adapter_Metrics.csv"),
        bclconvert_top_unknown = os.path.join(bclconvert_reports, "Top_Unknown_Barcodes.csv"),

        fastqc_read1 = expand(os.path.join(OUT_ROOT, "fastqc_run", "fastqc", "{samples}_R1_001_fastqc.zip"), samples = FIDs),
        fastqc_read2 = expand(os.path.join(OUT_ROOT, "fastqc_run", "fastqc", "{samples}_R2_001_fastqc.zip"), samples = FIDs),

        bbduk_log = expand(os.path.join(LOG_ROOT, "bbduk.{samples}.log"), samples = FIDs),

        fastqc_filtered_read1 = expand(os.path.join(OUT_ROOT, "fastqc_run", "fastqc", "{samples}_R1_bbduk_fastqc.zip"), samples = FIDs), 
        fastqc_filtered_read2 = expand(os.path.join(OUT_ROOT, "fastqc_run", "fastqc", "{samples}_R2_bbduk_fastqc.zip"), samples = FIDs), 

        bowtie2_R1 = expand(os.path.join(OUT_ROOT, "02_SILVA", "{samples}.DS.R1.bowtie2.log"), samples = FIDs),
        bowtie2_R2 = expand(os.path.join(OUT_ROOT, "02_SILVA", "{samples}.DS.R2.bowtie2.log"), samples = FIDs),

        kraken2_R1 = expand(os.path.join(OUT_ROOT, "03_kraken2", "{samples}.DS.R1.nt.report.kraken2"), samples = FIDs),
        kraken2_R2 = expand(os.path.join(OUT_ROOT, "03_kraken2", "{samples}.DS.R2.nt.report.kraken2"), samples = FIDs),

        bowtie2_genome_R1 = expand(os.path.join(OUT_ROOT, "02_REF", "{samples}.DS.genome_alignment.bowtie2.R1.log"), samples = FIDs),
        bowtie2_genome_R2 = expand(os.path.join(OUT_ROOT, "02_REF", "{samples}.DS.genome_alignment.bowtie2.R2.log"), samples = FIDs),

        multiQC_config = "resources/multiQC_config.yaml",
    output:
        multiQC = multiqc_report_path,
    conda:
        "multiqc"
    log:
        os.path.join(LOG_ROOT, "multiQC_report.log"),
    benchmark:
        os.path.join(BENCHMARK_ROOT, "multiQC_report.txt"),
    threads: 12
    resources:
        mem_gb = lambda wildcards, attempt: 24 + ((attempt - 1) * 2),
        time = lambda wildcards, attempt: 60 + ((attempt - 1) * 15),
    shell:
        "multiqc "
        "-n {multiqc_report_path} "
        "-s "
        "-f "
        "-c {input.multiQC_config} "
        "--interactive "
        "{input.bclconvert_runinfo} "
        "{input.bclconvert_demux_stats} "
        "{input.bclconvert_qual_metrics} "
        "{input.bclconvert_adpt_metrics} "
        "{input.bclconvert_top_unknown} "
        "{input.fastqc_read1} "
        "{input.fastqc_read2} "
        "{input.bbduk_log} "
        "{input.fastqc_filtered_read1} "
        "{input.fastqc_filtered_read2} "
        "{input.bowtie2_R1} "
        "{input.bowtie2_R2} "
        "{input.kraken2_R1} "
        "{input.kraken2_R2} "
        "{input.bowtie2_genome_R1} "
        "{input.bowtie2_genome_R2} "
        "2>&1 | tee {log}"



