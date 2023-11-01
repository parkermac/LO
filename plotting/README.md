# README for the plotting code

#### Flexible command line code for making snapshots or movies of model output.

---

#### `pan_plot.py`

This is the driver you use to make plots.  It can be run interactively or from the command line. It can make single snapshots, or movies over a time range.

---

#### `roms_plots.py`

A module of all the plot designs that are available through pan_plot. There is a hook in the code to look first for the existence of `LO_user/plotting/roms_plots.py` so you can easily customize your plots.

These make extensive use of the functions in the module `lo_tools/plotting_functions.py`.

---

#### `pinfo.py`

A module of scaling factors, colormaps, and other dicts that help keep the plot choices more centralized.

---

#### `create_sect_lines.p`

An interactive GUI tool for creating section lines, based on `LO/extract/tef2/create_sections.py`. At the command line you specify a grid to use (from your collection in LO_data/grids), and then you click to create one or more sections. These are saved in `LO_output/section_lines` as pickled DataFrames. The plotting function `roms_plots.P_sect_hc()` is the first instance of a plot that looks in that location for the section(s) it will use.
