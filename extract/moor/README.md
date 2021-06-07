# README for the mooring extractor code

## This code does mooring extractions and plotting.

---

`extract_moor.py` does the extractions. It uses `ncks` to get an extraction as a NetCDF file for each hour (or day), keeping them in a temporary folder and then concatenates them using `ncrcat`.  It is fast because it runs 100 ncks subprocesses at once. For command line arguments it relies on the shared tool `alpha/extract_argfun.py`.

`plot_moor.py` is a generic plotting tool for a user-selected mooring extraction.
