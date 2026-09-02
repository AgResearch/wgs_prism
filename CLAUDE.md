# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

`wgs_prism` is a Snakemake pipeline that turns a completed Illumina sequencing run into demultiplexed FASTQ plus a MultiQC report. It runs on AgResearch's **eRI** SLURM cluster against GPFS paths (`/projects`, `/datasets`, `/dataset`). It is operational infrastructure, not a library: there is no test suite, no CI, no linter config, and no package metadata.

The repository is small and every file in it is reachable from `run_wgs_qc`. A large body of unreachable code (a second workflow, six Python 2 scripts, a retired cluster profile, an unattended poller) was removed in the merge of PR #7 — see "What was removed" below before re-adding anything.

## Running it

Two layers of entry point, both shell:

| Script | Role |
|---|---|
| `run_wgs_qc` | Thin wrapper: resolves its own directory via `BASH_SOURCE`, `exec`s `_run_wgs_qc` — injecting `-i` only when called with no args, otherwise forwarding `"$@"` (injecting it unconditionally would make `./run_wgs_qc -r RUN` arrive as `-i -r RUN` and open the editor anyway) |
| `_run_wgs_qc` | The real driver (see below). `-i` interactive, `-r RUN` non-interactive (skips the editor) |

`_run_wgs_qc` does the work no Snakefile does: loads the modules, prompts for the run name, checks `RunInfo.xml` exists and infers the instrument type via `runinfo_to_machinetype.py` (erroring out if the type comes back empty), copies the run's sample sheet to `$OUT_ROOT/SampleSheet.csv`, and **when interactive opens it in `$EDITOR` (default `vim`) for hand-editing** (`OverrideCycles`, UMI settings, batch subsetting — the prompt text documents the cases), then invokes:

```sh
snakemake --profile config/slurm/eRI --snakefile workflow/wgs_prism.smk \
  --config OUT_ROOT=$OUTPUT_ROOT IN_ROOT=$SEQ_ROOT RUN=$RUN
```

Run that command directly to skip the interactive sample-sheet step. Add `-n` for a dry run — but note it cannot resolve the DAG past `checkpoint run_bclconvert` until real bclconvert output exists, so a dry run on a fresh run shows only the first rule.

Requires `module load Miniforge3/24.9.0-0` and `module load snakemake/7.32.3-foss-2023a-Python-3.11.6` first.

Both scripts resolve their own directory with `"$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"`, so a checkout runs its own code. `_run_wgs_qc` used to hardcode `/projects/2023_sequence_production/wgs_prism` and `cd` there, which made edits in any other checkout appear to have no effect.

There is no unattended poller in the repo. `-r` is the hook for one: it opens no editor and takes the sheet exactly as copied, so `_run_wgs_qc -r $RUN` is safe to drive from `cron`.

**Re-demultiplexing is blocked by design, in two independent places.** `_run_wgs_qc` exits if `$OUTPUT_ROOT/SampleSheet/bclconvert` already exists, and `run_bclconvert` refuses on its own: `--force` is deliberately **not** passed to bcl-convert. Because Snakemake pre-creates `bclconvert/` and `bclconvert/Logs/` from the declared outputs, the rule first `rmdir`s those empty shells; if real output is present the `rmdir` fails, the directory survives, and the rule aborts with an explicit message.

**The `rmdir` and the absent `--force` are one mechanism, not two independent choices.** Drop the `rmdir` while still omitting `--force` and bcl-convert fails on *every fresh run*, because Snakemake already created the directory. Add `--force` and the guard is gone. The remote cleanup branch went the other way (`--force`, no guard); the merge kept the guard deliberately. Change both or neither.

Because a retry can never succeed while that guard stands, `run_bclconvert` carries **`retries: 0`** and flat (non-`attempt`-scaled) `mem_gb=96`/`time=600`. Without it the profile's `restart-times` spent all five attempts on identical, doomed cluster jobs. Verified on a real sample-sheet error: `grep -c "Trying to restart job"` → 0.

Note what the guard does *not* mean. Snakemake removes the outputs of a failed job, and `bclconvert/` is a declared `directory()` output, so an ordinary rule failure **deletes the whole tree on the way out and a resubmit is not blocked** (observed: after a bcl-convert sample-sheet error the entire `SampleSheet/` directory was gone). Manual removal is needed in the cases where that cleanup never runs — the Snakemake controller itself killed or interrupted — or when you want to redo a demux that actually *succeeded*.

That same cleanup is why the rule captures `bcl_status=$?` instead of letting `set -e` kill the shell, then copies `bclconvert/Logs/*.log` to `$OUT_ROOT/logs/bclconvert_Logs/` and appends them to the rule log before `exit $bcl_status`. The logs of a failed demux live *inside* the directory Snakemake is about to delete, so that copy is the only surviving record. Keep the status capture — without it the copy is unreachable on exactly the path it was written for.

## Snakemake 7, not 8

`retries:` is a real rule keyword in 7.32.3 (`parser.py:446`) and overrides the profile's global `restart-times` per rule — `workflow.py:1642` falls back to the global only when `ruleinfo.retries is None`, so an explicit `retries: 0` genuinely disables restarts (`scheduler.py:738` tests `restart_times > attempt - 1`). `wgs_prism.smk` sets it on all three rules: `run_bclconvert: 0`, `run_fastqc: 3`, `run_multiqc: 1`.

**The profile no longer sets `restart-times` at all** (the cleanup removed it, and the only file that depended on the global default is gone). Every rule now carries an explicit `retries:`. A new rule that sets nothing gets one attempt — set it deliberately.

The eRI profile uses the **v7 cluster dialect**: `cluster:`, `cluster-status:`, `cluster-cancel:`, list-form `default-resources`. Applying Snakemake 8 knowledge (executor plugins, `--executor slurm`, `--jobs` semantics) will break `config/slurm/eRI/config.yaml`. `workflow/scripts/status.py` is the `cluster-status` hook — it maps `squeue`/`sacct` states to `success`/`failed`/`running`.

The profile's singularity bind mounts are what make the container-visible paths resolve:

```
/mnt/gpfs/persist/datasets:/datasets
/mnt/gpfs/persist/projects:/projects
/mnt/gpfs/scratch/projects:/scratch
/mnt/gpfs/persist/legacy_datasets:/dataset
```

So a `/datasets/...` reference path is a bind-mount path, not a host path.

## The workflow

**`workflow/wgs_prism.smk`** is the only workflow, and the only thing any entry point invokes. Three rules: `run_bclconvert` (a `checkpoint`) → `run_fastqc` per FASTQ → `run_multiqc` (declared `localrules`, so it runs on the submit host and the profile's cluster resources don't apply). Sample discovery is via `get_fastq_reports()`, which `os.listdir()`s the checkpoint's output directory at DAG-resolution time and filters out `Undetermined` — so samples come from files on disk, not from the sample sheet. bcl-convert is **not** a conda env — the rule does `export PATH=/agr/persist/apps/src/b/BCL-Convert:$PATH`. FastQC and MultiQC come from `workflow/envs/`.

`OUT_ROOT` names two different things and always has: `$OUT_ROOT/SampleSheet.csv` is the sample-sheet *file*, while `$OUT_ROOT/SampleSheet/` is the output *directory tree*. `wgs_prism.smk` uses `config["OUT_ROOT"]` raw and joins `"SampleSheet"` inline at each use site.

## Config, and the per-season edit

`config/pipeline_config.yaml` holds reference paths and defaults — after the cleanup, `multiqc_config` is the only reference key left. **`RUN`, `IN_ROOT` and `OUT_ROOT` are deliberately empty**, and the Snakefile refuses to start when any of them is missing *or* empty:

```python
for _required in ("RUN", "IN_ROOT", "OUT_ROOT"):
    if not config.get(_required):
        raise WorkflowError(f"{_required} must be set, e.g. --config {_required}=...")
```

They must be empty for that check to mean anything: the file used to ship a real-looking `RUN: "231005_A01439_0216_AHL7NKDRX3"` and a bare `OUT_ROOT`, so a forgotten `--config` key silently resolved against a stale run instead of stopping. Do not put example values back — the commented examples above each key are the record.

When the sequencing project rolls over (roughly annually), **one line changes**: `ILLUMINA_PROJECT` in `_run_wgs_qc` (~line 17), which drives `SEQ_ROOT` and `SEQ_BCLCONVERT_ROOT`. `IN_ROOT`/`OUT_ROOT` in `config/pipeline_config.yaml` are empty and no longer track the season; only the illustrative comments beside them mention a year.

Check `account:` in `config/slurm/eRI/config.yaml` at the same time. It is `2023_sequence_production`, which is live and deliberately **not** the same name as `ILLUMINA_PROJECT`; it does not follow the season.

Two reference paths remain outside the config file, which matters when updating either:

- `multiqc_config` → `config/pipeline_config.yaml`
- bcl-convert → `PATH` export inside the `run_bclconvert` rule body

`resources/multiQC_config.yaml` lists `bclconvert` and `fastqc` in `module_order` — exactly what the workflow produces — and carries `TBD` for the subtitle and sequencing platform.

## What was removed

The merge of PR #7 (`cleanup`) deleted ~4000 lines that nothing reached. `git show d3d3dac:<path>` recovers any of it. Do not reintroduce a file from this list without a reason:

- `workflow/detailed_follow_up.smk` — a deeper QC workflow (seqtk → bbduk → bowtie2/SILVA → kraken2/nt → bowtie2/GRCh38 → its own MultiQC). No entry point called it; its two `bowtie2_SILVA_alignment_*` rules were mangled by a snakefmt accident (string fragments dedented to column 0, falling out of the `shell:` expression); its GRCh38 index is gitignored, so a fresh clone could not run it anyway.
- `autostart_wgs_qc` — polled the retired `/dataset/2024_illumina_sequencing_d/active/` tree.
- `config/slurm/legacy/` — the profile for the retired HPC.
- Six Python 2-era `workflow/scripts/` files inherited from `gbs_prism` (`data_prism.py`, `kmer_prism.py`, `kmer_plots.r`, `reconcile_contaminants.py`, `add_sample_sheet_header.py`, `samplesheet_to_fastqname.py`), plus `resources/sample_sheet_header.csv` and `resources/adapters.fa`, which only they and the removed workflow read.
- Six `workflow/envs/*.yaml` for the removed workflow's tools.
- `kraken2_index` and `human_genome_index` from `config/pipeline_config.yaml`, and `restart-times` from the eRI profile.

`workflow/scripts/` now holds exactly two files, both wired in: `runinfo_to_machinetype.py` (called by `_run_wgs_qc` to map instrument prefix → `miseq`/`iseq`/`novaseq`) and `status.py` (the SLURM `cluster-status` hook).

## Conda environments

`workflow/envs/fastqc-0.12.1.yaml` and `multiqc-1.17.yaml` are fully-pinned `conda env export` dumps, and `multiqc-1.17.yaml` carries a machine-specific `prefix:` line. `workflow/envs/README.md` states the preferred approach is centrally-managed conda-store environments rather than letting Snakemake build these. Treat the YAMLs as a record of what was used, not as portable specs.
