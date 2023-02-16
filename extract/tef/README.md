# README for the tef code

### This is code to do TEF extractions and processing. The code also works on the segments between TEF sections, for example to make volume-integrated budgets of volume, salt and salinity variance.

#### WORKFLOW OVERVIEW

In order to get to where you can run `flux_salt_budget.py` you need to go through two separate workflows because these prepare the independent information at (1) TEF sections, and (2) segments between the sections.

(1) TEF sections:
 - `extract_sections.py`
 - `process_sections.py`
 - `bulk_calc.py`

(2) Segments (the flux_ code):
 - `make_segment_volumes.py`
 - `extract_segments.py`

Then you can run `tracer_budget.py`

**NOTE: this has not yet been updated to use xarray instead of netCDF4.**

---
#### PREAMBLE

`LO/pgrid/tef_section_maker.py` is a tool to define sections, if needed.  This is a GUI where you define section endpoints (forced to be either perfectly E-W or N-S), name the section, and define which direction is "landward" (i.e. which way is positive transport).  You can define one or more sections in a session.  At the end, when you push DONE it produces screen output suitable for pasting into new or existing lists in `tef_fun.get_sect_df()`.

`tef_fun.py` is an important module for this code. It includes `tef_fun.get_sect_df()` which returns a DataFrame of all the section names and their lon,lat endpoints.

---
#### EXTRACTION CODE

`extract_sections.py` creates a NetCDF file for each section with arrays of hourly transport and tracer values on the section, arranged as (t, z, x-or-y). Using command line arguments you can change the run, the day range, the sections to extract, and the variables extracted. Typically this will be done on a remote machine, like perigee, although the defaults are designed to work with model output I have saved on my mac.

**NOTE**: this code runs multiple subprocess instances of `extract_section_one_time.py`, (set by the Nproc command line argument which has a default value of 10). This significantly speeds things up, but it tends to occupy the machine, e.g. if you use -Nproc 20 on perigee you are using all the cores and may slow down other jobs.

**NOTE**: this code also automatically runs the two subsequent steps, `process_sections.py` and `bulk_calc.py`.  These can also be run as stand-alone (use -test True when running `extract_sections.py`) to facilitate debugging.

**PERFORMANCE**: For a full year of cas6, with -get_bio True and all NPZDOC variables this takes 5 hours on perigee (6.5 hours when including all steps).

The **command line arguments** are defined in `LO/lo_tools/lo_tools/extract_argfun.py`, with the usual requirements for gtagex, and beginning and end dates.  You also specify:
- -ro, --roms_out_num, is a integer specifying where to look for the ROMS history profiles
- -get_bio is a Boolean.  If True then you get the list in tef_fun.vn_list.  If False (the default) then vn_list = ['salt'].
- -sect_name is a string to specify the sections to get, either a single one such as ai1, or all of them using "-sect_name all" in the command line.  The full list is in tef_fun.get_sect_df().

Input: ROMS history files over some date range, e.g. [*] = 2017.01.01_2017.12.31

Output: LO_output/extract/[gtagex]/tef/extractions_[*]/[sect name].nc where:

Variables in the NetCDF files:
- salt is hourly salinity in each cell (t, z, x-or-y) [same for all other variables]
- q is hourly transport in each cell (t, z, x-or-y)
- vel is velocity in each cell (t, z, x-or-y) positive to East or North
- DA is the area of each cell (t, z, x-or-y) hence: q = vel * DA
- z0 is the average z-position of cell centers (assumes SSH=0), useful for plotting
- DA0 is the average cross-sectional area of each cell (assumes SSH=0)
- h is depth on the section (x-or-y) positive down
- zeta is SSH on the section (t, x-or-y) positive up
- ocean_time is a vector of time in seconds since (typically) 1/1/1970.

---

`process_sections.py` organizes all the transports at each time into salinity bins. It accepts a shortened set of command line arguments (-gatagex, -0, -1), or you can run it without these and it will ask you to choose from available folders.

Input: the output NetCDF files from `extract_sections.py`

Output: LO_output/extract/[gtagex]/tef/processed_[*]/[sect name].p

- a pickled dict with keys: ['q', 'salt', 'salt2', 'sbins', 'ot', 'qnet', 'fnet', 'ssh'] and other variables specified in tef_fun.vn_list

These are packed in order [time, salinity bin]:
- q is hourly transport in salinity bins (still positive East and North)
- salt is hourly transport of salt in salinity bins
- salt2 is hourly transport of salt-squared in salinity bins
- similar for all other variables that exist in both vn_list and the input file

The salinity bin centers (typically 1000 of them) are:
-	sbins

These remaining fields are just functions of time:
-	ot ocean time (sec from 1/1/1970)
-	qnet section integrated transport (m3/s)
-	fnet section integrated tidal energy flux (Watts)
-	ssh section averaged ssh (m)

---

`bulk_calc.py` does the TEF bulk calculations, using the algorithm of Marvin Lorenz, allowing for multiple in- and outflowing layers. Like `process_sections.py` it accepts a shortened set of command line arguments (-gatagex, -0, -1), or you can run it without these and it will ask you to choose from available folders.

Input: output of `process_sections.py`

Output: LO_output/extract/[gtagex]/tef/bulk_[*]/[sect name].p

These are pickled dicts with keys: ['salt', 'q', etc., 'ot', 'qnet', 'qabs' 'fnet', 'ssh'] and other variables, where 'salt', 'q', and etc. are matrices of shape (363, 30) meaning that it is one per day, at Noon UTC, after tidal-averaging, with nan-days on the ends cut off.  The 30 is the number of "bulk" bins, so many will be filled with nan's.

Since we have done the tidal averaging **qnet** is now the mean total transport through the section, mostly representing river flow in simple cases, whereas the qnet from `process_section.py` retains tidal variability. We also introduce a new time series, **qabs**, which is the tidal average of the absolute value of qnet.  The amplitude of the tidal flow is Qt = (pi/2)*qabs, and that Qprism = qabs/2 (from Chen et al. 2012 JPO).

---
#### TEF plotting

`bulk_plot.py` plots the results of bulk_calc.py, either as a single plot to the screen, or multiple plots to png's.  You need to edit the code to run for other years.

`check_results.py` is code with hard-wired file-paths to check that all out calculations above give exactly the same results as the original LiveOcean version. RESULT: the results are identical.

---
#### Code for the segments between TEF sections

`flux_fun.py` is a module of functions associated with the complexity of hooking up sections, segments and rivers for the salt budget and box model codes. Important functions are:
- flux_fun.get_two_layer() which sums the multi-layer output of `bulk_calc.py` into two layer TEF format and adjusts the signs to be positive in the direction of mean transport of the saltier layer. It is used in many of the plotting and analysis programs above.

---

`make_segment_volumes.py` gets the volume (with SSH = 0) of each segment.  Uses a cute mine-sweeper-like algorithm to find all the unmasked rho-grid cells on the grid between the sections that define a segment.

Input: flux_fun.segs, the model grid (hard wired in the code to be cas6), and sect_df = tef_fun.get_sect_df().

Output: LO_output/extract/tef/volumes_[gridname]/ containing:
- volumes.p a pickled DataFrame with segment names for the index and columns: ['volume m3', 'area m2', 'lon', 'lat']
- bathy_dict.p pickled dict of the grid bathy, for plotting (h, lon_psi, lat_psi)
- j_dict.p and i_dict.p pickled dict (keys = segment names) of the j and i indices of each segment area, stored as arrays ready for fancy indexing

---

`extract_segments.py` creates time series of hourly volume, salinity and other tracers averaged in segments, over a user specified (like a year) period.  Takes about 4 hours for a year on perigee.

Input: ROMS history files

Output: LO_output/extract/[gtagex]/tef/segments_[*].nc, an xarray Dataset whose variables have dimensions (time, segment).  Time is hourly.

---
#### BUDGET code

`tracer_budget.py` makes complete volume, and tracer budgets for a user-specified set of segments.

Input: TEF and segment-averaged tracer extractions, also some river extractions.

Output: Budget time series plots (one for each tracer) and some screen output of annual averages.
