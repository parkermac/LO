# README for LO/obsmod

### This code is designed to make comparisons between model runs and observations, for the purpose of model evaluation and to suggest areas of model improvement. It is meant to simplify this often-cumbersome task, e.g. making a full-year comparison of bottle observations from multiple sources into a single step. This code is based in earlier versions developed in LPM/obsmod.

---

**one_step_bottle_val_plot.py** This is a driver which acts on a model run [gtagex] for a chosen year or years. It does bottle cast extractions using code in LO/extract/cast for a series of sources (specified where?). It then combines the cast extractions with bottle data into a single pandas DataFrame, and then plots the result as a property-property plot, saving it as a png. Since it saved the figure to a file, it can be run directly on apogee, or whatever machine has direct access to model output files.

