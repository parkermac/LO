# README for the mooring extractor code

## This code does mooring extractions and plotting. It makes use of **xarray** instead of netCDF4 for manipulating things. The mooring extractions are packed in NetCDF files, with all 3-D data variables having dimensions ('ocean_time', 's_rho') or ('ocean_time', 's_w').

## Note that because of how xarray handles the time axis, it converts **ocean_time** to dtype='datetime64[ns]'.  See `plot_moor.py` for how to convert this to other things.

---

`extract_moor.py` does the extractions. It uses `ncks` to get an extraction as a NetCDF file for each hour (or day), keeping them in a temporary folder and then concatenates them using `ncrcat`.  It is fast because it runs 100 ncks subprocesses at once. For command line arguments it relies on the shared tool `alpha/extract_argfun.py`.

`plot_moor.py` is a generic plotting tool for a user-selected mooring extraction.


`pcolor_plot.py` is another generic plotting tool, that makes pcolor plots instead of lines.

`multi_mooring_driver.py` is a driver to run extract_moor.py for multiple moorings. It looks in Ldir['data']/moor at `job_lists.py` for a dict of station names and lon,lat tuples.  The idea of putting this in Ldir['data'] is that, while it is not in the LO repo, it is under the complete control of any user, meaning that the driver can be used by anyone directly, without the need for a user instance of anything. The downside is that you have to scp your own job_lists.py.
