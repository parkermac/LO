# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

This is the mooring extractor module of the **LiveOcean (LO)** ocean modeling framework. It extracts virtual mooring time-series from ROMS model output (NetCDF files) and plots the results. The broader LO project lives at `~/Documents/LO/` and its output at `~/Documents/LO_output/`.

## Running the scripts

**Quick test (single mooring):**
```bash
# In ipython:
run extract_moor -gtx cas7_t0_x4b -test True

# From command line:
python extract_moor.py -gtx cas7_t0_x4b -ro 0 -0 2017.07.04 -1 2017.07.06 -lt hourly -sn test -lon ' -125' -lat 47 -get_all True > test.log &
```
Note: negative longitude requires a leading space in quotes, e.g. `-lon ' -125'`.

**Multiple moorings at once:**
```bash
python multi_mooring_driver.py -gtx cas6_v3_lo8b -ro 2 -0 2019.07.04 -1 2019.07.06 -lt hourly -job orca -get_all True > mmd.log &
```

**Plotting:**
```bash
python plot_moor.py          # interactive file chooser, line plots
python plot_moor_vel.py      # same but focused on velocities (Godin lowpass applied)
python pcolor_plot.py        # pcolor plots of salt/temp/oxygen vs depth and time
```

## Key arguments

| Flag | Meaning |
|------|---------|
| `-gtx` | `gridname_tag_exname`, e.g. `cas7_t0_x4b` |
| `-ro` | `roms_out_num` (0 = default roms_out) |
| `-0`/`-1` | start/end date as `YYYY.MM.DD` |
| `-lt` | `hourly`, `daily`, `lowpass`, `average` |
| `-sn` | station name (used in output filename) |
| `-Nproc` | parallel `ncks` subprocesses (default 10; mac/perigee handle 100 fine) |
| `-get_tsa/vel/bio/surfbot/pressure` | variable group flags |
| `-get_all` | shorthand for all of tsa+vel+bio+surfbot |

## Architecture

### Extraction pipeline (`extract_moor.py`)
1. Parses the ROMS grid to find the nearest wet grid cell (`find_good()`) — checks rho, u, v grids separately to handle masked land cells.
2. Spawns many `ncks` subprocesses in batches of `Nproc` to extract one file per time step into a temp folder.
3. Concatenates with `ncrcat` into a single output NetCDF at `LO_output/extract/{gtagex}/moor/{sn}_{ds0}_{ds1}.nc`.
4. Adds computed `z_rho`/`z_w` vertical coordinate arrays via xarray, then copies attributes from the source files.

### Multi-mooring driver (`multi_mooring_driver.py`)
Calls `extract_moor.py` sequentially for each station in a `sta_dict`, then moves results into `LO_output/extract/{gtagex}/moor/{job_name}/`. Logs per-station stdout/stderr to `…/moor/logs/`.

### Station lists (`job_lists.py`)
`get_sta_dict(job_name)` returns `{station_name: (lon, lat)}`. The driver looks for this file first in `LO_user/extract/moor/` and falls back to the copy here in `LO/extract/moor/`. Add new jobs as `elif` branches in `get_sta_dict`.

### Output format
Output NetCDF files use xarray conventions. `ocean_time` is stored as `datetime64[ns]`. 3-D variables have dimensions `('ocean_time', 's_rho')` or `('ocean_time', 's_w')`; 2-D variables have just `('ocean_time',)`. Converting time for plotting:
```python
ot_dt = pd.to_datetime(ds.ocean_time.values)
t = (ot_dt - ot_dt[0]).total_seconds().to_numpy()
```

### Dependencies
- `lo_tools` (LO framework utilities: `Lfun`, `zrfun`, `zfun`, `plotting_functions`)
- NCO tools: `ncks` (extraction), `ncrcat` (concatenation) — must be on PATH
- Python: `xarray`, `numpy`, `matplotlib`, `pandas`, `cmocean`
