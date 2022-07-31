# README for LO/obs/dfo

#### This code is for working with the Canadian CTD and bottle data that was provided by Elise Olson, a former postdoc in Susan Allen's group at UBC. More complete notes about the original files are in LO_data/obs/dfo.

---

`process_ctd_bottle.py` is the main code to run to do the complete processing. The input/output naming follows the general LO/obs format.

---

`loadDFO.py` is a module of functions for reading in selected parts of the sqlite databases in LO_data/obs/dfo.

`test_loadDFO.py` is some code to test the parts of `loadDFO.py`, and to start simple explorations of the data.

`test_process_[ctd,bottle].py` is code used to further develop the and test the two parts of the processing.

`plot_[ctd,bottle].py` plots the processing results for a selected year.
