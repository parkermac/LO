# README for LO/extract/tide

## This code is used only to quickly extract hourly tide height time series from a model run. The same job could be done with extract_moor.py, but that by default gets the full water column and works on one station location at a time. This should be much faster.

#### extract_tide.py

This created pickled Series of model tide height at stations where we have observations.

Input: station numbers and locations from LO_output/obs/tide/info_[station_number]_[year].csv

Output: pickled Series in LO_output/extract/[gtagex]/tide/tide_[station_number]_[year].p

Algorithm strategy: Get grid indices for all stations, then for each model history file get the zeta field and extract all the locations at once using fancy indexing. Save all Series to individual pickle files at the end.
