# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

`wgs_prism` is a Snakemake pipeline that turns a completed Illumina sequencing run into demultiplexed FASTQ plus a MultiQC report. It runs on AgResearch's **eRI** SLURM cluster against GPFS paths (`/projects`, `/datasets`, `/dataset`). It is operational infrastructure, not a library: there is no test suite, no CI, no linter config, and no package metadata.

## Running it

Three layers of entry point, all shell:

| Script | Role |
|---|---|
| `run_wgs_qc` | Thin wrapper: resolves its own directory via `BASH_SOURCE`, `exec`s `_run_wgs_qc` — injecting `-i` only when called with no args, otherwise forwarding `"$@"` (injecting it unconditionally would make `./run_wgs_qc -r RUN` arrive as `-i -r RUN` and open the editor anyway) |
| `_run_wgs_qc` | The real driver (see below). `-i` interactive, `-r RUN` non-interactive (skips the editor) |
| `autostart_wgs_qc` | Polls every 15 min for `RTAComplete.txt`, then calls `_run_wgs_qc -r $RUN` |

`_run_wgs_qc` does the work no Snakefile does: loads the modules, prompts for the run name, infers the instrument type via `runinfo_to_machinetype.py`, copies the run's sample sheet to `$OUT_ROOT/SampleSheet.csv`, **opens it in `$EDITOR` (default `vim`) for hand-editing** (`OverrideCycles`, UMI settings, batch subsetting — the prompt text documents the cases), then invokes:

```sh
snakemake --profile config/slurm/eRI --snakefile workflow/wgs_prism.smk \
  --config OUT_ROOT=$OUTPUT_ROOT IN_ROOT=$SEQ_ROOT RUN=$RUN
```

Run that command directly to skip the interactive sample-sheet step. Add `-n` for a dry run — but note it cannot resolve the DAG past `checkpoint run_bclconvert` until real bclconvert output exists, so a dry run on a fresh run shows only the first rule.

Requires `module load Miniforge3/24.9.0-0` and `module load snakemake/7.32.3-foss-2023a-Python-3.11.6` first.

**Re-demultiplexing is blocked by design, in two independent places.** `_run_wgs_qc` exits if `$OUTPUT_ROOT/SampleSheet/bclconvert` already exists, and `run_bclconvert` refuses on its own: `--force` is deliberately **not** passed to bcl-convert. Because Snakemake pre-creates `bclconvert/` and `bclconvert/Logs/` from the declared outputs, the rule first `rmdir`s those empty shells; if real output is present the `rmdir` fails, the directory survives, and the rule aborts with an explicit message. Bypassing the wrapper does not get you a rerun of a *successful* demux — remove `$OUT_ROOT/SampleSheet/bclconvert` by hand for that.

Because a retry can never succeed while that guard stands, `run_bclconvert` carries **`retries: 0`** and flat (non-`attempt`-scaled) `mem_gb=96`/`time=600`. Without it the profile's `restart-times: 5` spent all five attempts on identical, doomed cluster jobs. Verified on a real sample-sheet error: `grep -c "Trying to restart job"` → 0.

Note what the guard does *not* mean. Snakemake removes the outputs of a failed job, and `bclconvert/` is a declared `directory()` output, so an ordinary rule failure **deletes the whole tree on the way out and a resubmit is not blocked** (observed: after a bcl-convert sample-sheet error the entire `SampleSheet/` directory was gone). Manual removal is needed in the cases where that cleanup never runs — the Snakemake controller itself killed or interrupted — or when you want to redo a demux that actually *succeeded*.

## Snakemake 7, not 8

`retries:` is a real rule keyword in 7.32.3 (`parser.py:446`) and overrides the profile's global `restart-times` per rule — `workflow.py:1642` falls back to the global only when `ruleinfo.retries is None`, so an explicit `retries: 0` genuinely disables restarts (`scheduler.py:738` tests `restart_times > attempt - 1`). `wgs_prism.smk` sets it on all three rules: `run_bclconvert: 0`, `run_fastqc: 3`, `run_multiqc: 1`. The profile's `restart-times: 5` now only supplies the default for `detailed_follow_up.smk`, whose rules set nothing.

The eRI profile uses the **v7 cluster dialect**: `cluster:`, `cluster-status:`, `cluster-cancel:`, list-form `default-resources`. Applying Snakemake 8 knowledge (executor plugins, `--executor slurm`, `--jobs` semantics) will break `config/slurm/eRI/config.yaml`. `workflow/scripts/status.py` is the `cluster-status` hook — it maps `squeue`/`sacct` states to `success`/`failed`/`running`.

The profile's singularity bind mounts are what make the container-visible paths resolve:

```
/mnt/gpfs/persist/datasets:/datasets
/mnt/gpfs/persist/projects:/projects
/mnt/gpfs/scratch/projects:/scratch
/mnt/gpfs/persist/legacy_datasets:/dataset
```

So `kraken2_index: /datasets/...` is a bind-mount path, not a host path.

`config/slurm/legacy/` is the profile for the retired HPC (`account=perrybe`, `inv-*` partitions, `restart-times: 0`, `cluster-status` commented out). Nothing uses it; edit `config/slurm/eRI/` unless you specifically mean the old cluster.

## The two workflows

**`workflow/wgs_prism.smk`** — the production pipeline, and the only one any entry point invokes. Three rules: `run_bclconvert` (a `checkpoint`) → `run_fastqc` per FASTQ → `run_multiqc` (declared `localrules`, so it runs on the submit host and the profile's cluster resources don't apply). Sample discovery is via `get_fastq_reports()`, which `os.listdir()`s the checkpoint's output directory at DAG-resolution time and filters out `Undetermined`. bcl-convert is **not** a conda env — the rule does `export PATH=/agr/persist/apps/src/b/BCL-Convert:$PATH`.

**`workflow/detailed_follow_up.smk`** — 511 lines of deeper QC (seqtk downsample → bbduk trim → bowtie2 vs SILVA → kraken2 vs nt → bowtie2 vs GRCh38 → its own MultiQC). **Not reachable from any entry point.** Two things to know before touching it:

- It discovers samples with `glob_wildcards()` on `bclconvert/{sample}_R1_001.fastq.gz`, evaluated at **parse time**. It must be run *after* `wgs_prism.smk` has produced FASTQ, or `FIDs` is empty and the DAG silently does nothing.
- The `shell:` blocks of `bowtie2_SILVA_alignment_read1`/`_read2` are mangled: the `"-U ..."`, `"1> /dev/null"`, `"2> {output}"` string fragments are dedented to column 0, so they fall out of the `shell:` expression and become no-op module-level statements. The command that actually runs is `bowtie2 -p N -x SILVA138.1` with no input and no output redirect — so the rule fails on missing output rather than producing a bad log. Looks like a snakefmt/black accident. Fix before relying on these rules.

Its target rule is `default`; invoke with `--snakefile workflow/detailed_follow_up.smk` and the same `--config` triple.

**`OUT_ROOT` means two different things in the two Snakefiles.** `wgs_prism.smk` uses `config["OUT_ROOT"]` raw and joins `"SampleSheet"` inline at each use site; `detailed_follow_up.smk` rebinds a module-level `OUT_ROOT = os.path.join(config["OUT_ROOT"], "SampleSheet")` near the top. The resolved paths agree, but copying a path expression between the two files will double or drop the `SampleSheet` segment. Compounding it: `$OUT_ROOT/SampleSheet.csv` is the sample-sheet *file*, while `$OUT_ROOT/SampleSheet/` is the output *directory tree*.

## Config, and the per-season edit

`config/pipeline_config.yaml` holds reference paths and defaults. **`RUN`, `IN_ROOT` and `OUT_ROOT` are deliberately empty**, and both Snakefiles refuse to start when any of them is missing *or* empty:

```python
for _required in ("RUN", "IN_ROOT", "OUT_ROOT"):
    if not config.get(_required):
        raise WorkflowError(f"{_required} must be set, e.g. --config {_required}=...")
```

They must be empty for that check to mean anything: the file used to ship a real-looking `RUN: "231005_A01439_0216_AHL7NKDRX3"` and a bare `OUT_ROOT`, so a forgotten `--config` key silently resolved against a stale run instead of stopping. Do not put example values back — the commented examples above each key are the record.

When the sequencing project rolls over (roughly annually — see the last two commits), **one line changes**: `ILLUMINA_PROJECT` in `_run_wgs_qc` (~line 17), which drives `SEQ_ROOT` and `SEQ_BCLCONVERT_ROOT`. `IN_ROOT`/`OUT_ROOT` in `config/pipeline_config.yaml` are empty and no longer track the season; only the illustrative comments beside them mention a year.

Check `account:` in `config/slurm/eRI/config.yaml` at the same time. It is `2023_sequence_production`, which is live and deliberately **not** the same name as `ILLUMINA_PROJECT`; it does not follow the season.

Reference-data paths are inconsistently located, which matters when updating an index:

- `human_genome_index`, `kraken2_index` → `config/pipeline_config.yaml`
- SILVA 138.1 → **hardcoded** in `detailed_follow_up.smk` (carries a literal `# TODO parameterize in config`)
- bcl-convert → `PATH` export inside the `run_bclconvert` rule body

`resources/GRCh38/` is gitignored, so a fresh clone has no human genome index and `detailed_follow_up.smk`'s alignment check cannot run.

## Live vs. inherited code

`workflow/scripts/` is mostly **dead inheritance from the `gbs_prism` pipeline**. Only two files are wired in:

- `runinfo_to_machinetype.py` — called by `_run_wgs_qc` to map instrument prefix → `miseq`/`iseq`/`novaseq`
- `status.py` — the SLURM `cluster-status` hook in the eRI profile

The other six (`data_prism.py`, `kmer_prism.py`, `kmer_plots.r`, `reconcile_contaminants.py`, `add_sample_sheet_header.py`, `samplesheet_to_fastqname.py`) are referenced by **nothing** in this repo — as is `resources/sample_sheet_header.csv`, whose only consumer is `add_sample_sheet_header.py -H`. They are Python 2-era (`from __future__ import print_function`, bare `reduce`, `sys.path.append('/dataset/gseq_processing/active/bin/gquery')`) and hardcode legacy `/dataset/gseq_processing/...` paths. Don't read them to understand the pipeline, and don't assume they run.

`autostart_wgs_qc` is likewise stale: its `get_run_folder` probes `/dataset/2024_illumina_sequencing_d/active/` three times identically, while `_run_wgs_qc` now reads from `/projects/$ILLUMINA_PROJECT/run_data`. It will not find current runs without editing. Its *second* blocker is now fixed: `_run_wgs_qc` used to open the editor unconditionally, so even a run it did find would hang forever on `vim` with no terminal. The editor is now inside `if [ "$INTERACTIVE" = yes ]`, and `-r` proceeds with the sheet exactly as copied.

## Conda environments

`workflow/envs/*.yaml` are fully-pinned `conda env export` dumps, and some carry a machine-specific `prefix:` line (e.g. `multiqc-1.17.yaml` points at `/home/.../perrybe/.conda/envs/multiqc`). `workflow/envs/README.md` states the preferred approach is centrally-managed conda-store environments rather than letting Snakemake build these. Treat the YAMLs as a record of what was used, not as portable specs.
