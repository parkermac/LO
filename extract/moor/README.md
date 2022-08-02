# README for the mooring extractor code

## This code does mooring extractions and plotting.

It makes use of xarray instead of netCDF4 for manipulating things. The mooring extractions are packed in NetCDF files, with all 3-D data variables having dimensions ('ocean_time', 's_rho') or ('ocean_time', 's_w'), or just 'ocean_time'.

Note that because of how xarray handles the time axis, it converts **ocean_time** to dtype='datetime64[ns]'.  See `plot_moor.py` for how to convert this to other things.

---

#### Extraction

`extract_moor.py` does the extractions. It uses `ncks` to get an extraction as a NetCDF file for each hour (or day), keeping them in a temporary folder and then concatenates them using `ncrcat`.  It is fast because it runs many ncks subprocesses at once (use the optional -Nproc flag to control this).  My mac and perigee can handle 100 fine although it may slow other operations.  Boiler seemed to choke on 100 but ran fine with 10.

`multi_mooring_driver.py` is a driver to run extract_moor.py for multiple moorings. It looks in `LO_user/extract/moor/job_lists.py` for a dict of station names and (lon,lat) tuples. If that `job_lists.py` file does not exist, then it uses the one in this directory.

See the codes for details on the required command line arguments. Basically you need to tell it which run to use, and the time limits and frequency. For `extract_moor.py` you also pass a station name, longitude, and latitude, whereas for `multi_mooring_driver.py` you instead pass a job name.

---

#### Plotting

`plot_moor.py` is a generic plotting tool for a user-selected mooring extraction.

`pcolor_plot.py` is another generic plotting tool, that makes pcolor plots instead of lines.

`plot_moor_standalone.py` is a standalone version of `plot_moor.py` that I give to people when sending them mooring extractions, to help them get started.
