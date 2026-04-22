# README for pcRaft

Author: Dakota Mascarenas

Last updated: 2026/04/13

These files are Penn Cove mussel raft sonde data from the work conducted and discussed in, or adjacent to, Roberts & Carrington (2023): https://doi.org/10.1016/j.jembe.2023.151927

Data received via email on 2026/01/02 and 2026/02/24 from Emily Carrington to Dakota Mascarenas, with metadata clarification also on 2026/02/24.

Location: 48.220861°N, 122.705667°W

Two stations are processed:
* `raft_main`: Primary YSI sonde data from 2014-2019 (hourly). Variables: temperature, salinity, chlorophyll, pH, DO. Source file: `PennCove-mussel-raft-sonde-data-2014-2019.xlsx` (sheet: `cleaned_data`).
* `raft_hiRes`: July-December 2017 high-resolution HOBO temperature sondes at 0.5-7 m depth (originally 10-minute resolution, resampled to hourly). Source files: `B8_[depth]m.csv`. Temperature only — salinity is not directly measured.

NOTE: Data type is `moor` (mooring format) despite being measured by YSI sondes.

NOTE: Original timestamps are in PST (UTC-8); converted to UTC during processing.

NOTE: `raft_hiRes` does not have directly measured salinity associated with temperature values. CT is computed using salinity interpolated from `raft_main` at the nearest time and depth using the GSW toolbox.

NOTE: GSW conversions applied are from practical salinity (SP) to absolute salinity (SA) and in-situ temperature to conservative temperature (CT) using the `gsw` library.

NOTE: DO converted from mg/L to uM (x1000/32). Chl are in ug/L = mg/m3.

Mooring data availability (`raft_main`):
* CT: 2014-2019
* SA: 2014-2019
* DO: 2014-2019
* PH: 2014-2019
* Chl: 2014-2019

Mooring data availability (`raft_hiRes`):
* CT: 2017