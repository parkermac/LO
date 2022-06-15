# README for the plotting code

#### Flexible command line code for making snapshots or movies of model output.

---

#### `pan_plot.py`

This is the driver you use to make plots.  It can be run interactively or from the command line.

#### `roms_plots.py`

A module of useful plot designs.  These make extensive use of the functions in the module `lo_tools/plotting_functions.py`. There is now a hook in the code to look first for the existence of `LO_user/plotting/roms_plots.py`.

#### `pinfo.py`

A module of scaling factors, colormaps, and other dicts that help keep the plot choices more centralized.
