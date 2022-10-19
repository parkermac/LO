# README for hycom

## This code is designed to make a backfill archive for the GOFS 3.0 and 3.1 versions of HYCOM.

NOTE: HYCOM has a 22-year reanalysis 1994-2015: GLBv0.08, 53x.  I have not made use of this yet but it could be great for longer hindcasts.  The notes below are for more recent years, and are used for the current system.

---

#### A typical workflow would consist of two steps.  First, from the command line do:

```
python get_dt_lists.py
```

If you were doing it for the first time when building an archive from scratch you would do
`python get_dt_lists.py -a True` (which forces it go get everything, WARNING: This. Takes. Days.).

Next, actually get any required hycom extractions by doing:

```
python get_hycom_days.py > log_days &
```

Use -test True to test, and -force to ensure new files, e.g. after testing. Otherwise the default is only to get the new files that are not already there.

---

`hfun.py` is a module of functions.

The list of grids and experiments we choose from is in hfun.hy_dict.  For example, GLBy0.08 is a grid with 0.08 deg lon x 0.04 deg lat that covers 80S to 90N.  The older GLBu0.08 is 2x coarser in latitude.

The geographic area for extractions is set in hfun.get_extraction_limits().

---

`get_dt_lists.py` creates lists of available daily datetimes for each experiment, saving them in LO_data/hycom/dt_lists/hy1.p, for example.  Associated info is in dt_info_hy1.txt.

The default is that this just updates the list for the last experiment in hy_dict.  You can update all by running it with the flag -a True.

**As mentioned above, you should run this first when updating the archive of extractions, because the extraction code relies on the lists to specify which days to look for.**

---

`get_hycom_days.py` it the main driver for filling out the archive.  While it can get daily extractions for the full set of experiments in hfun.hy_dict, in practice it is much more parsimonious, only getting an extraction if it is not already in the archive.  It works on all days in the dt_lists.

Output: is the raw NetCDF file that comes from HYCOM (NetCDF4).  These include the usual coordinate and metadata, so they are the preferred format.  In later processing we account for the needed transformations to what ROMS expects.  This is done for example in ocn0.  More notes on details are in test_extraction.py.

The results are stored, for example, in: **LO_data/hycom/hy1/h2012.01.25.nc**, one for each entry in the associated dt_list.

Also outputs a file log.txt with more compact info on success or failure.

---

`print_dt_info.py` is a convenience function to print all of the dt_info files in one place.  Here is typical output:

```
Experiment = hy1
 - glb = GLBu0.08, exnum = 90.9
 - 2012.01.25 to 2013.08.20
 - missing 133 out of 573 days

Experiment = hy2
 - glb = GLBu0.08, exnum = 91.0
 - 2013.08.17 to 2014.04.08
 - missing 0 out of 234 days

Experiment = hy3
 - glb = GLBu0.08, exnum = 91.1
 - 2014.04.07 to 2016.04.18
 - missing 4 out of 742 days

Experiment = hy4
 - glb = GLBu0.08, exnum = 91.2
 - 2016.04.18 to 2018.11.20
 - missing 26 out of 946 days

Experiment = hy5
 - glb = GLBu0.08, exnum = 93.0
 - 2018.09.20 to 2018.12.09
 - missing 5 out of 80 days

Experiment = hy6
 - glb = GLBy0.08, exnum = 93.0
 - 2018.12.05 to 2021.04.25
 - missing 14 out of 872 days
```

---

`test_extraction.py` is code I used for testing different methods of getting a single backfill day.
