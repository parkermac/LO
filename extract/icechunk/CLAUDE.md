# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Environment

This code runs in the `loenv` conda environment (defined at `LO/loenv.yml`). Key packages: `icechunk`, `virtualizarr`, `obstore`, `obspec_utils`, `xarray`, `zarr`.

Activate with: `conda activate loenv`

## Scripts

### `liveocean_icechunk_example.py` — virtual Icechunk store from forecast layers

```bash
# Dry run — print derived S3 URLs without reading or writing
python liveocean_icechunk_example.py --dry-run

# Create a new Icechunk store (default: 7 days before the known URL date)
python liveocean_icechunk_example.py --target liveocean-source --mode create

# Append to an existing store
python liveocean_icechunk_example.py --target liveocean-source --mode append --days-before 1

# Target source.coop instead
python liveocean_icechunk_example.py --target source-coop --mode create
```

### `rechunker_main.py` — SLURM job launcher for history file rechunking

Runs on klone (UW HPC). Launches `rechunker_worker.sh` as a SLURM array job, one task per timestep.

```bash
# Testing on mac
python rechunker_main.py -gtx cas7_t2_x11b -0 2024.01.01 -1 2024.01.02 -lt hourly0 -test True

# Production on klone
python rechunker_main.py -gtx cas7_t2_x11b -0 2024.01.01 -1 2024.01.01 -lt hourly0
python rechunker_main.py -gtx cas7_t1_x11ab -0 2024.01.01 -1 2024.01.31 -lt average
```

`list_type` choices: `hourly0` (first day: hours 0–24), `hourly` (hours 1–24), `average` (daily averages). Use `hourly0` when populating a new collection; switch to `hourly` when appending.

## Required environment variables

Credentials are loaded via `load_dotenv()` from `~/dotenv/<target>.env` (selected automatically by `--target`).

| Variable | Used for |
|---|---|
| `AWS_ACCESS_KEY_ID` | Write credentials for target bucket |
| `AWS_SECRET_ACCESS_KEY` | Write credentials for target bucket |
| `LIVEOCEAN_ICECHUNK_PREFIX` | S3 key prefix for `liveocean-source` target |
| `SOURCE_COOP_BUCKET` / `SOURCE_COOP_PREFIX` | For `source-coop` target |

The source bucket (`liveocean-share` on kopah) requires signed credentials — kopah disables anonymous access at the daemon level even for public buckets.

## Architecture

### Virtual Icechunk pipeline (`liveocean_icechunk_example.py`)

No data is copied. VirtualiZarr references point back to the original `.nc` files on S3.

**Data flow:**
1. `build_source_urls()` — derives `layers.nc` S3 paths for N days before a known reference date using the `f{YYYY.MM.DD}` folder convention.
2. `build_virtual_dataset()` — calls `open_virtual_dataset()` (VirtualiZarr + HDFParser) per file, loading only `ocean_time` into memory. Standardizes time coords (`ocean_time` → `time` + `step`) and concatenates along `time`.
3. `write_icechunk()` — creates or opens an Icechunk `Repository` on S3, writes virtual references via `ds.virtualize.to_icechunk()`, and commits.

**Storage targets** are configured in `get_target_config()`: `liveocean-source` (kopah, path-style S3) or `source-coop` (AWS-style). The `VirtualChunkContainer` in `make_icechunk_storage()` tells Icechunk how to resolve chunk reads back to the kopah source.

### Rechunker pipeline (`rechunker_main.py` + `rechunker_worker.sh` + `rechunker_worker_one_time.py`)

Designed for klone/SLURM. Converts ROMS history files from kopah S3 into rechunked NetCDF files in `/gscratch`, then (later) ingests into Icechunk.

**Data flow:**
1. `rechunker_main.py` — builds date ranges, submits SLURM array jobs via `sbatch --array=1-N rechunker_worker.sh`.
2. `rechunker_worker.sh` — SLURM batch script; activates `loenv`, sets `$TMPDIR` per task, calls `rechunker_worker_one_time.py`.
3. `rechunker_worker_one_time.py` — maps `$SLURM_ARRAY_TASK_ID` to a specific history file, copies it to `$TMPDIR` via `s5cmd`, subsets variables with `ncks`, writes to `/gscratch` output directory.

Source files are read from `s3://liveocean-<username>/LO_roms/<gtagex>/`. Output goes to `$LOo/icechunk_temp/sub_<date>_<list_type>/`.
