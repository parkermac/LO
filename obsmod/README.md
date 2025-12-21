# README for LO/obsmod

### This code is designed to make comparisons between model runs and observations, for the purpose of model evaluation and to suggest areas of model improvement. It is meant to simplify this often-cumbersome task, e.g. making a full-year comparison of bottle observations from multiple sources into a single step. This code is based in earlier versions developed in LPM/obsmod.

---

**one_step_bottle_val_plot.py** This is a driver which acts on a model run [gtagex] for a chosen year or years. It does bottle cast extractions using **LO/extract/cast/extract_casts_fast.py** for a list of sources. It then combines the cast extractions with bottle data into a single pandas DataFrame, and then plots the result as a property-property plot, saving it as a png. Since it saved the figure to a file, it can be run directly on apogee, or whatever machine has direct access to model output files.

It relies on these two stand-alone programs:

**combine_obs_mod.py** Code to combine observed and modeled bottle values for a collection
of sources.

**plot_val.py** Makes a figure with a bunch of property-property plots and a map, comparing a model run with observations. Note that if you download a pickled DataFrame created by combine_obs_mod.py then you can make use of several optional command line arguments to make more focused plots, e.g. plotting only casts inside the Salish Sea.

**obsmod_functions.py** is a module of things that might be useful for more than one piece of code here. This is where the list of observational data sources to loop over is specified.

