# wgs_prism

Demultiplex and QC a completed Illumina sequencing run.

The pipeline reads a run directory from the sequencer and gives one MultiQC report. It
demultiplexes with Illumina `bcl-convert`. MGI DNBSEQ runs are different. Do not use this
pipeline for them. Use [`mgi_prism`](../mgi_prism).

The pipeline does these tasks:

1. It demultiplexes the run with `bcl-convert` and the sample sheet that you edit.
2. It runs FastQC on each demultiplexed fastq file.
3. It collects the `bcl-convert` reports and the FastQC reports into one MultiQC report.

## Before you start

You need these three items:

* **The run.** The run directory must be under `/projects/$ILLUMINA_PROJECT/run_data/<RUN>`. It
  must contain `RunInfo.xml`. The launcher reads the instrument from this file and selects the
  machine type. The run is complete when `RTAComplete.txt` is present.
* **A sample sheet.** The launcher uses `<RUN>/SampleSheet.csv`. If that file is absent, it
  accepts one `*.csv` file in the run directory whose name occurs in the run name. Copy the sheet
  into the run directory before you start.
* **An output directory.** The default is `/projects/$ILLUMINA_PROJECT/postprocessing/illumina/<machine
  type>/<RUN>`. The machine type is `novaseq`, `miseq` or `iseq`. Give a different directory for a
  second batch of the same run.

## Operation

Start the launcher from the pipeline directory:

```sh
./run_wgs_qc
```

The launcher has eight steps. It submits no job before step 8.

| Step | Task | Notes |
|---|---|---|
| 1 | Give the run name | For example `260831_A01439_0491_AHVWLGDRX7`. The launcher confirms that the directory is present. |
| 2 | Confirm the output directory | The default is `<POSTPROCESSING>/<machine type>/<RUN>`. The launcher makes the directory. |
| 3 | Confirm the run and the output directory | The launcher shows both paths. |
| 4 | Confirm a run that is not complete | The launcher asks this question only when `RTAComplete.txt` is absent. The default answer is **no**. |
| 5 | The launcher copies the sample sheet | It copies to `OUT_ROOT/SampleSheet.csv`. It copies only when that file is absent. |
| 6 | Edit the copy | The launcher opens `$EDITOR`, or `vim` when `$EDITOR` is not set. See "How to edit the sample sheet". |
| 7 | Confirm the edited sample sheet | The launcher then confirms that `OUT_ROOT/SampleSheet/bclconvert` is absent. |
| 8 | Confirm | The launcher shows the command, then submits to SLURM. |

The launcher changes the source sample sheet at no time. It edits the copy in the output
directory.

Step 5 keeps a copy that is already there. A second launch into the same output directory
therefore gives you the sheet that you edited before. This behaviour helps you continue a run
that failed before the demultiplex. It does not give you a second batch, because step 7 stops a
launch into a directory that holds demultiplexed output. Delete `OUT_ROOT/SampleSheet.csv` to
start again from the sheet in the run directory.

The launcher writes the output of snakemake to `OUT_ROOT/logs/wgs_prism.log`. It shows no
progress on the screen. Read that file to see the state of the run.

### Batches

A batch is the set of rows that you keep in the sample sheet at step 6. A batch can have any
size, down to one sample. Use a batch when one run holds libraries that need different
`OverrideCycles` settings.

**To demultiplex a second batch of the same run, start the launcher again and give a different
output directory**, for example `.../260831_A01439_0491_AHVWLGDRX7_BATCH2`. The launcher then
copies the sheet from the run directory again, and you keep the rows of the second batch. A
second launch into the same output directory cannot operate, because that directory holds the
demultiplexed output of the first batch.

The `bcl-convert` reports give the percentages of all the clusters on the flowcell. A small batch
therefore shows a large `Undetermined` fraction in `Demultiplex_Stats.csv` and in the
`bcl-convert` section of the MultiQC report. This result is correct, and it is not an error of
the batch.

**A run that holds more than one index length must go into more than one batch.** `bcl-convert`
accepts one `OverrideCycles` setting for each run, and that setting gives one index length. Count
the index lengths before you start:

```sh
awk '/^\[Data\]/{d=1;next} d&&/,/{print}' SampleSheet.csv | \
  awk -F, 'NR==1{for(i=1;i<=NF;i++)if($i=="index")c=i;next}{print length($c)}' | sort | uniq -c
```

More than one line in that result means more than one batch. Give each batch its own output
directory and its own `OverrideCycles` line. A sheet that mixes the lengths stops with:

```
ERROR: Sample Sheet Error: SampleSheet.csv sample #97 (index 'GTCAAGTCCA') has an index of
length 10 bases, but a length of 8 was expected based upon RunInfo.xml or the OverrideCycles
setting.
```

### How to edit the sample sheet

At step 6, do these tasks:

* **Delete the rows that are not part of this batch.** A second batch of the same run needs a
  second launch and a different output directory.
* **Add an `OverrideCycles` line under `[Settings]`** when the batch needs one. The table gives
  the usual cases.

| Case | `[Settings]` lines |
|---|---|
| 151 bp PE, 10 bp i5, no i7 | `OverrideCycles,Y151;I10;N10;Y151` |
| 101 bp PE, 8 bp indices | `OverrideCycles,Y101;I8N2;I8N2;Y101` |
| 151 bp PE, 8 bp indices | `OverrideCycles,Y151;I8N2;I8N2;Y151` |
| 101 bp SE, 10 bp i5, 8 bp i7 | `OverrideCycles,Y101;I10;I8N2` |
| T-overhang removal | `OverrideCycles,N1Y150;I10;I10;N1Y150` |
| 101 bp PE, 12 bp UMI, 8 bp i5 and i7 | `OverrideCycles,Y101;U12I8;I8;Y101` and `CreateFastqForIndexReads,1` and `TrimUMI,0` |
| 101 bp PE, 20 bp UMI, 8 bp i7 | `OverrideCycles,Y101;U20;I8;Y101` and `CreateFastqForIndexReads,1` and `TrimUMI,0` |

For a run that has UMIs: when the sheet gives the indices with `N` bases for the UMI, remove
those bases and keep the index sequence only.

References:
[T-overhang trimming](https://support.illumina.com/bulletins/2020/06/trimming-t-overhang-options-for-the-illumina-rna-library-prep-wo.html),
[UMI settings](https://knowledge.illumina.com/software/on-premises-software/software-on-premises-software-reference_material-list/000007337).

## Warnings and what to do

| Message | Meaning | Action |
|---|---|---|
| `sorry can't find the sample-sheet for this run` | The run directory has no `SampleSheet.csv` and no `*.csv` file whose name occurs in the run name. | Copy the sheet into the run directory. The launcher stops. |
| `RTAComplete.txt does not exist` | The run is still in progress, or still in transfer. | Wait for the transfer to end. The default answer is no. |
| `Warning: <directory> already exists, use anyway ?` | The output directory holds the results of an earlier run. | Snakemake keeps the complete results and makes the remainder. Use this behaviour to continue a run that failed. |
| `<OUT_ROOT>/SampleSheet/bclconvert already exists` | The run is demultiplexed. | See "Re-demultiplexing" below. The launcher submits nothing. |
| `already contains demultiplexed output; refusing to overwrite` | The same condition, from the workflow. | See "Re-demultiplexing" below. |
| `has an index of length N bases, but a length of M was expected` | The sample sheet holds more than one index length, or `OverrideCycles` does not agree with the indices. | Count the index lengths, then give each length its own batch. See "Batches". |

### Re-demultiplexing

Two independent guards refuse to write over demultiplexed output. The launcher stops when
`OUT_ROOT/SampleSheet/bclconvert` is present, and the `run_bclconvert` rule stops for the same
reason. The pipeline does not give `--force` to `bcl-convert`.

**A retry of the demultiplex cannot succeed while that output is there.** The rule therefore
gives `retries: 0`, and it makes one attempt only. Its memory and its time limit are single
values that are large enough for the worst case, and they do not increase at each attempt.
Before this, the pipeline made five identical attempts that each failed in the same way.

**An ordinary failure does not block the next attempt.** Snakemake removes the output of a rule
that fails, and `bclconvert` is a declared output directory. A demultiplex that stops on a bad
sample sheet, or on a job that the scheduler kills, therefore removes its own output. Correct
the sample sheet and start the launcher again.

Remove the directory by hand in two conditions only: the demultiplex was successful and you want
to do it again, or snakemake itself stopped before it removed the output.

```sh
rm -rf $OUT_ROOT/SampleSheet/bclconvert
```

## Outputs

The pipeline writes these files into the output directory:

| Path | Contents |
|---|---|
| `SampleSheet/multiqc/<RUN>.multiqc.html` | The report. Open this file first. |
| `SampleSheet/multiqc/<RUN>.multiqc_data/` | The tables of the report, in TSV format. |
| `SampleSheet/bclconvert/<sample>_S<n>_L00<n>_R<n>_001.fastq.gz` | The demultiplexed reads, and `Undetermined`. |
| `SampleSheet/bclconvert/Reports/` | The reports of `bcl-convert`. `Demultiplex_Stats.csv`, `Quality_Metrics.csv` and `Top_Unknown_Barcodes.csv` are the usual ones. |
| `SampleSheet/bclconvert/Logs/` | The logs of `bcl-convert`, and `FastqComplete.txt`. |
| `SampleSheet/fastqc_run/fastqc/` | The FastQC report of each fastq file. |
| `SampleSheet.csv` | The copy of the sample sheet that defines this batch. |
| `logs/wgs_prism.log` | The full output of snakemake. |
| `logs/run_bclconvert.log`, `logs/run_fastqc.<fastq>.log` | The log of each step. |
| `logs/slurm/` | The output of each SLURM job. |
| `benchmarks/` | The run time and the memory of each job. |

The MultiQC report gives one FastQC entry for each fastq file. A paired-end run of four lanes
thus gives eight entries for each sample. The pipeline excludes the `Undetermined` reads from
FastQC. The `bcl-convert` section of the report includes them.

## Options

The file `config/pipeline_config.yaml` holds the reference paths.

| Key | Default | Function |
|---|---|---|
| `multiqc_config` | `resources/multiQC_config.yaml` | The MultiQC configuration. It sets the title, the logo and the sequence of the sections. |
| `kraken2_index` | `/datasets/2024-kraken2-indices/k2_nt_20231129` | The Kraken2 index. `detailed_follow_up.smk` uses it only. |
| `human_genome_index` | `/projects/2023_sequence_production/wgs_prism/resources/GRCh38/GRCh38` | The bowtie2 index of GRCh38. `detailed_follow_up.smk` uses it only. |

**The keys `RUN`, `IN_ROOT` and `OUT_ROOT` in that file are empty.** Each workflow refuses to
start when one of them has no value. The launcher gives a value for each of them with `--config`.
Keep them empty: a value that looks real lets a forgotten `--config` key read, or write over, a
different run.

The `kraken2_index` path is a bind-mount path in the container, and not a path on the host. The
profile makes these paths available. See "Environment" below.

`resources/multiQC_config.yaml` gives `TBD` for the subtitle and for the sequencing platform.
Change these two values when you want them in the report. Its `module_order` also lists bbduk,
SILVA, kraken and GRCh38. `wgs_prism.smk` makes no input for those modules, and MultiQC omits
those sections without a message. `detailed_follow_up.smk` makes them.

## Unattended operation

Load the modules first. The version of snakemake is a fixed requirement. The profile
`config/slurm/eRI` uses version 7 syntax, and later versions reject it.

```sh
module purge && module load Miniforge3/24.9.0-0
module load snakemake/7.32.3-foss-2023a-Python-3.11.6
snakemake --profile config/slurm/eRI --snakefile workflow/wgs_prism.smk \
  --config OUT_ROOT=$OUT_ROOT IN_ROOT=$SEQ_ROOT RUN=$RUN
```

This command does not open an editor. Use it when the sample sheet at
`$OUT_ROOT/SampleSheet.csv` is already correct.

**Give all three keys.** The workflow stops immediately when one of them is absent or empty:

```
WorkflowError in file workflow/wgs_prism.smk, line 21:
RUN must be set, e.g. --config RUN=...
```

`detailed_follow_up.smk` makes the same check. Earlier versions held a value for each key in
`config/pipeline_config.yaml`, and a forgotten key then resolved against that value without a
message.

**`snakemake -n` gives an incomplete answer on a new run.** The rule `run_bclconvert` is a
checkpoint. Snakemake finds the samples only after that rule writes its output. A dry run on a
run that is not demultiplexed thus shows the first rule only.

---

# For developers

`CLAUDE.md` gives the structure between the files and the invariants of the domain.

## Layers

| File | Function |
|---|---|
| `run_wgs_qc` | A short wrapper. It finds its own directory and runs `_run_wgs_qc`. It gives `-i` when you give it no option, and it gives your options to `_run_wgs_qc` when you give one. |
| `_run_wgs_qc` | The interactive launcher. It loads the modules, selects the run and the output directory, copies the sample sheet, opens the editor, and runs snakemake. |
| `workflow/wgs_prism.smk` | The workflow. Three rules. |
| `workflow/scripts/runinfo_to_machinetype.py` | It maps the instrument prefix to the machine type. `M` gives `miseq`, `F` gives `iseq`, and each other letter gives `novaseq`. |
| `workflow/scripts/status.py` | The `cluster-status` hook of the profile. It maps the states of `squeue` and `sacct` to `success`, `failed` or `running`. |

**Each script finds its own directory**, with
`WGS_PRISM_BIN="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"`. A checkout thus runs its own
code. `_run_wgs_qc` held the literal `/projects/2023_sequence_production/wgs_prism` before, and a
checkout in a different directory therefore ran the production code, and an edit appeared to have
no result. `mgi_prism` and `gtseq_prism` use the same method.

**The editor is inside `if [ "$INTERACTIVE" = yes ]`.** The flag `-r` gives the sheet to
`bcl-convert` exactly as copied, and it opens no editor.

`autostart_wgs_qc` polls for `RTAComplete.txt` and then calls `_run_wgs_qc -r $RUN`. **It does not
operate at present.** Its `get_run_folder` probes `/dataset/2024_illumina_sequencing_d/active/`
three times, and that path is retired. Correct that path before you use it. Its second condition
is now correct: the editor no longer stops the unattended path.

## Workflow structure

The rules run in this sequence: `run_bclconvert`, then `run_fastqc` for each fastq file, then
`run_multiqc`.

`run_bclconvert` is a **checkpoint**. The function `get_fastq_reports()` lists the output
directory of that checkpoint when snakemake resolves the DAG, and it removes the `Undetermined`
files. Snakemake thus finds the samples from the files on disk, and not from the sample sheet.

`run_fastqc` makes one job for each fastq file. Read 1 and read 2 are two jobs.

`run_multiqc` is in `localrules`. It runs on the submit host. The resources of the profile do not
apply to it.

Each rule gives its own `retries` value. This keyword is available in snakemake 7.32.3, and it
replaces the `restart-times` of the profile for that rule.

| Rule | `retries` | Reason |
|---|---|---|
| `run_bclconvert` | 0 | A second attempt finds the output of the first one and stops in the same way. One attempt only, with a memory and a time limit that are large enough for the worst case. |
| `run_fastqc` | 3 | This rule gives the same result at each attempt, and its failures are usually transient. The memory and the time limit increase at each attempt. |
| `run_multiqc` | 1 | The last step, after some hours of work. A failure here is usually permanent, and one more attempt is sufficient. |

The `restart-times: 5` of the profile now applies to `detailed_follow_up.smk` only, because the
rules of that file give no `retries` value.

`bcl-convert` is not a conda environment. The rule adds `/agr/persist/apps/src/b/BCL-Convert` to
`PATH`. FastQC and MultiQC come from `workflow/envs/`.

### The two guards on re-demultiplexing

Snakemake makes `bclconvert/` and `bclconvert/Logs/` before the rule runs, because the rule
declares them as outputs. The rule therefore runs `rmdir` on these two empty directories first.
The `rmdir` fails when the directories hold real output, the directory survives, and the rule
stops with a message. Keep this sequence when you change the rule.

## Environment

* **Snakemake 7.32.3 is a fixed requirement.** The profile `config/slurm/eRI` uses version 7
  syntax (`cluster:`, `cluster-status:`, `cluster-cancel:` and a list-form `default-resources`).
  Version 8 replaces this syntax with executor plugins and refuses the profile.
* **The SLURM account is in the profile.** The file `config/slurm/eRI/config.yaml` sets
  `account: 2023_sequence_production` in `default-resources`. This value is the current
  allocation. It is independent of `ILLUMINA_PROJECT`, and the two names are not the same.
* **The bind mounts of the profile make the container paths resolve.** The profile binds
  `/mnt/gpfs/persist/datasets` to `/datasets`, `/mnt/gpfs/persist/projects` to `/projects`,
  `/mnt/gpfs/scratch/projects` to `/scratch`, and `/mnt/gpfs/persist/legacy_datasets` to
  `/dataset`.
* `config/slurm/legacy/` is the profile of the retired HPC. Nothing uses it.
* The files `workflow/envs/*.yaml` are `conda env export` dumps, and some hold a `prefix:` line of
  one machine. `workflow/envs/README.md` gives the preferred approach: central conda-store
  environments. Treat these files as a record and not as portable recipes.

## The per-season edit

The sequencing project changes each year. **One line changes:** `ILLUMINA_PROJECT` in
`_run_wgs_qc`. It sets `SEQ_ROOT` and `SEQ_BCLCONVERT_ROOT`. The keys `IN_ROOT` and `OUT_ROOT` in
`config/pipeline_config.yaml` are empty, and they do not follow the season.

Confirm `account:` in `config/slurm/eRI/config.yaml` at the same time. It does not follow
`ILLUMINA_PROJECT`, and the two names are not the same.

## detailed_follow_up.smk

`workflow/detailed_follow_up.smk` is a deeper QC workflow: seqtk downsample, bbduk trim, bowtie2
against SILVA 138.1, kraken2 against nt, bowtie2 against GRCh38, and its own MultiQC report. **No
entry point calls it.** Its target rule is `default`.

Three conditions apply before you use it:

* **Run it after `wgs_prism.smk`.** It finds the samples with `glob_wildcards()` when snakemake
  reads the file. An empty `bclconvert/` directory thus gives an empty DAG and no error.
* **The two `bowtie2_SILVA_alignment_*` rules are damaged.** The fragments `"-U ..."`,
  `"1> /dev/null"` and `"2> {output}"` are at column 0. They are outside the `shell:` expression,
  and they operate as statements of the module. The command that runs has no input and no output
  redirect. Correct these rules first.
* **The SILVA index is a literal** at `/datasets/2024-silva-rrna/SILVA138.1`. The file carries a
  `# TODO parameterize in config` comment. The directory `resources/GRCh38/` is in `.gitignore`,
  and a new clone therefore has no human genome index.

**`OUT_ROOT` means two different things in the two workflow files.** `wgs_prism.smk` uses
`config["OUT_ROOT"]` and adds `"SampleSheet"` at each use. `detailed_follow_up.smk` sets a module
variable `OUT_ROOT = os.path.join(config["OUT_ROOT"], "SampleSheet")`. The paths agree, but a path
expression that moves between the two files gains or loses the `SampleSheet` part. Also,
`OUT_ROOT/SampleSheet.csv` is the sample sheet **file**, and `OUT_ROOT/SampleSheet/` is the output
**directory**.

## Inherited code

Six files in `workflow/scripts/` come from the `gbs_prism` pipeline: `data_prism.py`,
`kmer_prism.py`, `kmer_plots.r`, `reconcile_contaminants.py`, `add_sample_sheet_header.py` and
`samplesheet_to_fastqname.py`. The file `resources/sample_sheet_header.csv` is with them. **No
file in this repository calls them.** They are Python 2 code, and they hold paths of the retired
`/dataset/gseq_processing/` tree. Do not read them to understand the pipeline.

There is no test suite, no CI and no linter configuration. This repository is operational
infrastructure.
