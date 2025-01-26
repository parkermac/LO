# README for LO/obs/tide

## This is code to download and process tide height data.

### This is the starting point for testing tide height simulation in the model. The whole process consists of three parts:
- Get observed tide time series (this code)
- Extract the corresponding time series from a model run. This is handled in LO/extract/obs, and relies on the sta_df files created here in order to know the station locations.
- Comparisons are handled in LPM/obsmod/compare_tide.py.

## There are two programs to run here. You have to edit them by hand to get different years.

### get_noaa_tide.py
Input: from a NOAA website

Output: LO_output/obs/tide/
- sn_df_noaa_[year].p a pickled DataFrame of station names and locations. The index is the unique station number.
- sn_df_noaa_[year].csv a csv version for convenience of inspection
- tide_noaa_[station number]_[year].p pickled Series of a year of SSH vs. time, hourly, UTC, at a given station.

### get_dfo_tide.py
Input: from a DFO website

Output: LO_output/obs/tide/
- sn_df_dfo_[year].p a pickled DataFrame of station names and locations. The index is the unique station number.
- sn_df_dfo_[year].csv a csv version for convenience of inspection
- tide_dfo_[station number]_[year].p pickled Series of a year of SSH vs. time, hourly, UTC, at a given station.
