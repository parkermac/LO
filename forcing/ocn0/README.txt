These programs create the ocean files used to force a roms run, typically for a single day or a 3 day forecast.  They bring in hycom data, either archived extractions for backfill, or from a web source for a forecast.

It has been substantially recoded from previous versions to use the new hycom FMRC_best format of the forecast files.

For backfill it will use an archive of files of similar format, available from late 2012 to the present, created by the code in hycom2.

2019.05.21 The new method Ofun.get_interpolated_alt() makes extensive use of nearest neighbor interpolation and the results is that the code is substantially faster: 1 min vs. 7 min for a cas6 backfill day.

2020.03.16 I added some new methods to Ofun.py to make it get HYCOM FMRC files one day at a time.

2020.08.13 Added a method to get the HYCOM files using the nco operator "ncks"

More update notes in the separate pieces of code.

===============================================================================
* make_forcing_main.py is the main driver, similar in basic construction and usage to all the code of the same name in other forcing folders (e.g. riv2, tide1, atm).

Input: command line arguments tell it what to do - specifically which gridname_tag = [gtag] to work on, what day [f_string], and whether this is a forecast or backfill.  Then it either gets data from the web, or from LiveOcean_data/hycom2.

Output: LiveOcean_output/[gtag]/[f_string]/ocn4/...
(1) Info/process_status.csv and screen_out.txt (great for debugging)
(2) Data/...
- h[data_string].nc (all the hycom data for this period - only for the new FMRC_best forecast)
- coord_dict.p (pickled dict of coordinate vectors - just like in hycom1)
- h[data_string].p (pickled dict of fields for a day - for a forecast there will be 8 of these, 6 for backfill)
- fh[date_string].p (time filterd version of h[].p files - for a forecast there will be 4 of these, 2 for backfill)
- xfh[date_string].p (spatially extrapolated version of fh[].p files - for a forecast there will be 4 of these, no nans)
(3) ocean_[clm, ini, bry].nc are the files for forcing roms

Notes:
- If something goes wrong with getting the data for a forecast, then it uses "planB" in which the ocean_clm.nc file is copied from the day before, and one day is added to its last time.
- in backfill mode it uses the archives files , e.g., LiveOcean_data/hycom2/hy6/h2019.01.01.nc.  Currently no planB for this operation.

Modules:
- Ofun.py: the main workhorse functions to get data, filter in time, and extrapolate
- Ofun_bio.py: add and fill bio variables, e.g. using regressions against salinity
- Ofun_CTD.py: fill fields in the Salish Sea and coastal estuaries using CTD observations.  Currently only set to work for 1/1/2017.  The goal is to get a better initial condition.
- Ofun_nc.py: functions for making the ouput NetCDF files (3).

===============================================================================
* check_ocn.py plot ocean fields to compare raw hycom, extrapolated hycom, and interpolated roms.

Input: a selection of the files for a given [f_string]/ocn4, as well as the grid.nc file from LiveOcean_data/grids.

Output: a screen plot of map fields for a variable - useful for debugging

===============================================================================
* check_bio.py plot ocean fields to check on interpolated roms fields of the bio variables.

Input: ocean_ini.nc a given [f_string]/ocn4, as well as the grid.nc file from LiveOcean_data/grids.

Output: a screen plot of map fields for a bunch of variables - useful for debugging

===============================================================================
* compare_ocn4_ocn4old.py makes comparisons of the results of ocn4old and ocn4, for a specific day and grid.  The only difference between the two is that ocn4 uses Ofun.get_interpolated_alt() whereas ocn4old uses Ofun.get_interpolated(), and this difference is implemented around line 184 of make_forcing_main.py.

Input: ocean_ini.nc for a specific [gtag]/[f_string]/ocn4old and /ocn4, as well as the grid.nc file from LiveOcean_data/grids.

Output: a series of screen plots of map fields for all variables, comparing the two calculations.  You can change whether it plots the surface or bottom layer using a switch in the code.

RESULT: the comparisons looked close enough to indicate that the new fast interpolation is working as planned, and so I have gotten rid of ocn4old.
