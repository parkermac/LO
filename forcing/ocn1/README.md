# README for the ocn1 forcing

## This is just like ocn0 except I fixed a bug in the time units.

These programs create the ocean files used to force a roms run, typically for a single day or a 3 day forecast.  They bring in hycom data from archived extractions for backfill, or from a web source for a forecast.

For backfill it uses an archive of files of similar format, available from late 2012 to the present, created by the code in LO/pre/hycom and placed in LO_data/hycom.

---

## make_forcing_main.py
This is the main driver, similar in basic construction and usage to all the code of the same name in other forcing folders (e.g. riv0).

### Input:
Command line arguments tell it what to do - specifically which gridname_tag = [gtag] to work on, what day [f_string], and whether this is a forecast or backfill.  Then it either gets data from the web, or from LO_data/hycom.

### Output: (*) = LO_output/forcing/[gtag]/[f_string]/ocn0
- (*)/ocean_[clm, ini, bry].nc are the files for forcing roms
- (*)/Info/results.txt and screen_output.txt (great for debugging)
- (*)/Data/...
  - h[data_string].nc (all the hycom data for this period, one per day)
  - coord_dict.p (pickled dict of coordinate vectors)
  - h[data_string].p (pickled dict of fields for a day - for a three-day forecast there will be 8 of these, for a backfill day there will be 6)
  - fh[date_string].p (time filterd version of h[].p files - for a forecast there will be 4 of these, 2 for backfill)
  - xfh[date_string].p (spatially extrapolated version of fh[].p files - for a forecast there will be 4 of these, no nans)

### Notes:
- If something goes wrong with getting the data for a forecast (planA and planB), then it uses "planC" in which the ocean_clm.nc file is copied from the day before, and one day is added to its last time.
- in backfill mode it uses the archives files , e.g., LO_data/hycom/hy6/h2019.01.01.nc.  No planB for this operation.

### Modules:
- Ofun.py: the main workhorse functions to get data, filter in time, and extrapolate
- Ofun_bio.py: add and fill bio variables, e.g. using regressions against salinity
- Ofun_CTD.py: fill fields in the Salish Sea and coastal estuaries using CTD observations.  Currently only set to work for a specific day.  The goal is to get a better initial condition.
- Ofun_nc.py: functions for making the ouput NetCDF files.

---

## check_results.py
- still need to add this.
