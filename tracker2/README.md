# README for tracker

## particle tracking utility for the LO system

This code is designed to be a flexible, and hopefully fast, tool for doing particle tracking experiments using saved history files from ROMS.  It assumes that you save fields hourly, and that the file structure follows LiveOcean standards (hours 1-24, numbered 00[02-25]) in a single folder named by date, e.g. f2019.12.18.  This code uses nearest neighbor for most everything and so might work with more complex ROMS grids (the LO grids are plaid), but that is untested.

---

#### Steps to run a particle tracking experiment:

(1) Create the nearest neighbor search trees for your grid.  You only need to do this once for a grid, but likely will want to repeat on different machines or after updating python.  You can do this from the command line with something like:
```
python make_KDTrees.py -gtx cas7_t0_x4b -d 2017.07.04 -ro 0
```
but providing the tags appropriate for your run.  This takes a few minutes and creates some pickled KDTrees in LO_output/tracker2_trees/[gridname]/.

---

(2) Edit `LO_user/tracker2/experiments.py` to define an experiment.

To create the initial positions of your release you make new entries in:

- get_ic() where you define the lon, lat, and fractional depth vectors for the initial conditions associated with an experiment name.  The examples like "jdf0" are the simplest ones to copy from as they just involve making a regular lon, lat grid of points near the surface. Particles on land are trimmed automatically by the tracker code.

---

(3) Run `tracker.py` from the command line, with arguments that tell it which experiment to use, the start times, and some other choices like whether or not to track in 3-D.  It should be run from the linux command line because, for reasons I don't understand, the performance gets slower after repeated runs in ipython.  A typical run command might be like:
```
python tracker.py -exp jdf0 -3d True -d 2019.7.04 -dtt 2
```
which would track in 3d (including turbulence by default) for two days, starting on 2019.07.04 at a bunch of locations near the surface in the Strait of Juan de Fuca.  Important choices are encoded in the name of the output directory.  If it finds an existing directory of the same name it will warn you and quit.  If you want to overwrite that directory, add "-clb True" (clobber).

Look at the code near the top of `tracker.py` to see all possible arguments and their default values.  In a single experiment, for example, you can have many start days, separated by any number of days.  You can also add an optional tag to the end of the experiment folder name.

The output appears as NetCDF files in, for example:

LO_output/tracks2/[gtagex]/jdf0_3d/
- exp_info.csv (a record of the experiment choices)
- grid.nc (bathy info, used for plotting)
- release_2019.07.04.nc (there would be more than one of these if you used more start days)

NOTE: exp_info.csv is ALSO written to LO_output/tracks2 because we need some grid info when reading in things at the start of trackfun.py, i.e. when it is imported. This is NOT good coding because it can cause errors when trying to run two jobs at the same time. I never figured out how to pass arguments to modules, and maybe that would also be a bad idea.

NOTE: the default is to save output every hour, even though the underlying calculation may be done in finer steps (the default ndiv=12 means that we use 300 s steps in the RK4 integration).  You can save more frequently using the "sph" (saves per hour) input parameter.  I find this to be convenient during debugging.

---

(4) Plot the results using `tplot.py` as a first generic tool.  Copy this to a new name in LO_user/tracker2 to start making your own plotting tool.

Run from ipython on your laptop, because it creates a screen plot.  It will prompt you to choose the output you want to plot.

NOTE: In the code I set a parameter npmax=300 which subsamples the full set of tracks so there are no more than npmax.  This keeps the plotting from being too slow.

---

#### Under the hood there are some things you should not need to edit:

`trackfun_nc.py` handles creating and appending to the NetCDF output files

`trackfun.py` is the real heart of the program.  It is a module of functions that does all steps of an experiment, typically in one-day chunks as orchestrated by the calling code `tracker.py`.

**NOTE: `tracker.py` will automatically look first for "LO_user/tracker2/trackfun.py" which is a hook in case you want to make your own edited version of the functions to do something exotic like adding diurnal cycling behavior to particles.**

LIMITATIONS: Currently the code is hardwired to only save time series of particle positions, velocity, salinity and temperature.  I will need to do a bit more work to simplify adding more tracers.

---

#### New Development Notes (most recent on top)

2024.08.12:
- This version in the tracker2 folder has several bug fixes
compared to the original in the tracker folder. These became apparent when we
started doing tracking in very high-resolution runs in Hood Canal.
- We fixed the logic for updating the vertical position. It had been using the relative
depth of the starting location. Now instead we figure out how much z should change and use it to update relative depth at the end location for each time step.
- We added the ability (should be a flag argument) to use more nearest neighbors
to find the vertical velocity. Spikey w applied over 5 minute time steps caused
excessive vertical motion.
- We changed the nearest neighbor trees to use meters in all directions. For
high resolution grids the original lon, lat, pcs trees were unreliable.

2022.11.16:
- Based on experiments by Jilian Xiong we made two changes to the turbulent mixing in `trackfun.py`. (i) We modified the top and bottom AKs values to be the same as those one grid cell down (from the top) or up (from the bottom), so that the nearest neighbor did not get a near-zero AKs when it was near the boundaries. (ii) We modified the calculation of d(AKs)/dz to use the instantaneous (tidally varying) dz. Experiments showed that these made little difference, but performed slightly better in the "well-mixed" tests.

2022.10.06:
- added a new module called `customizations.py` that allows more complex behavior, such as switching to a different target depth after a certain number of days.

2021.10.14:
- added "tracer_list_full" near the top of `trackfun.py` to allow more tracers to be written to the output (and potentially used along the way).
- entirely done in xarray
- TO DO: `trackfun_nc.py` needs a way to automate adding tracer units and long names.  Or I could just add them all by hand.

2021.10.05:
- started refactoring for the LO code structure

---

#### Old Development notes from when this was in the LiveOcean repo (most recent on top)

2021.05.15:
- Added a new flag "sh" for start hour, allowing the release to begin on an arbitrary GMT hour of the first day.

2021.03.03:
- Got rid of old tracker folder and renamed this one tracker (was tracker2 in)

2020.11.10:
- Dropped the "_ndiv" from the output file name if the default (ndiv=12) is used.
- Added optional "sink" flag, e.g. use "-sink 40" to add 40 meters per day sinking speed to all particles.  Then the tag "_sink40" is added to the output file.  It would be nice to generalize to allow floating and sinking with less confusing naming (e.g. floating would be _sink-40 which is not intuitive).

2020.10.04:
- Introduced 3-point Hanning window (0.25, 0.5, 0.25) vertical smoothing on AKs, which improved WMC vmix results.  This fix was based on many experiments using test_dAKs_dz.py and modeled after work in North et al. (2006).

2020.10.03 and the past week:
- Updated experiments.py to simplify naming logic: drop ic_name, use exp_name instead.  Also cleaned up and added some new methods like ic_from_list() and ic_from_TEFsegs().
- Allow "cs" limits to go over full -1 to 0 range.  This works even for finding tracer values because of the nearest neighbor search.
- Changed the method for calculation of dAKs/dz guided by results of WMC (Well Mixed Condition) in vmix experiment.
- Introduced a new optional input flag for tracker.py.  "-no_advection True" turns off all advection, meaning only mixing is active.  Used with vmix experiment for WMC test.
- Also added the "-tag" input flag, allowing you to append a tag to an experiment folder.
