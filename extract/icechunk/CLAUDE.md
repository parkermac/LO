# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Environment

This code runs in the `loenv` conda environment (defined at `LO/loenv.yml`). Key packages: `icechunk`, `virtualizarr`, `obstore`, `obspec_utils`, `xarray`, `zarr`.

Activate with: `conda activate loenv`

## Running the script

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

## Required environment variables

Credentials are loaded from env files under `~/dotenv/`. The script calls `load_dotenv()` automatically (currently commented out — `from dotenv import load_dotenv` must be re-enabled if using `.env` files).

| Variable | Used for |
|---|---|
| `AWS_ACCESS_KEY_ID` | Write credentials for target bucket |
| `AWS_SECRET_ACCESS_KEY` | Write credentials for target bucket |
| `LIVE_OCEAN_TARGET` | Set automatically from `--target` arg |
| `LIVEOCEAN_ICECHUNK_PREFIX` | S3 key prefix for `liveocean-source` target |
| `SOURCE_COOP_BUCKET` / `SOURCE_COOP_PREFIX` | For `source-coop` target |

The source bucket (`liveocean-share` on kopah) is accessed anonymously via `skip_signature=True` / `anonymous=True`. The target bucket always requires credentials.

## Architecture

`liveocean_icechunk_example.py` builds a **virtual** Icechunk store: no data is copied. It creates VirtualiZarr references that point back to the original `.nc` files on S3, then commits those references into an Icechunk repository.

**Data flow:**
1. `build_source_urls()` — derives `layers.nc` S3 paths for N days before a known reference date using the `f{YYYY.MM.DD}` folder convention.
2. `build_virtual_dataset()` — calls `open_virtual_dataset()` (VirtualiZarr + HDFParser) for each file, loading only `ocean_time` into memory. Standardizes time coordinates (`ocean_time` → `time` + `step`) and concatenates along `time`.
3. `write_icechunk()` — creates or opens an Icechunk `Repository` on S3, writes virtual references via `ds.virtualize.to_icechunk()`, and commits with a descriptive message.

**Storage targets** are configured in `get_target_config()`: `liveocean-source` (kopah, path-style) or `source-coop` (AWS-style). The `VirtualChunkContainer` inside `make_icechunk_storage()` tells Icechunk how to resolve chunk reads back to the kopah source.

## Known issues / caveats

- `TARGET_STORAGE_NAME` is referenced in `parse_args()` but commented out at the top — the script will raise `NameError` unless `--storage-name` is passed or the constant is restored.
- `load_dotenv` is imported but commented out; the `main()` call to it will fail unless re-enabled.
- Kopah (`s3.kopah.uw.edu`) requires signed credentials even for public buckets — anonymous access is disabled at the daemon level by UW sysadmins.
