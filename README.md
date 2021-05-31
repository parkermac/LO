# README for LO

## LO is the code base for running the LiveOcean collection of regional ocean simulations.

The code here handles all the non-research tasks of making model (e.g. ROMS) forcing files, running the model on a remote linux cluster, and post-processing of model output.  Many of the post-processing tasks, such as mooring- cast- and layer extractions, and particle tracking, are designed for other users to be able to run from the linux command line, or modify if needed.

This version is the second major release of this project, the first was called LiveOcean. This version incorporates all the things I learned building the first one. A few specifics: it uses python Path objects everywhere possible. It defaults to using NetCDF4 for all the ROMS forcing. It uses python code instead of shell scripts wherever that makes things simpler - especially in the drivers. It also uses a much more rigorous and naming and organization system.

All the instructions assume you are working from the linux (bash) command line.

Contact Parker MacCready, p.maccready@gmail.com, with any questions.

Check out today's model output at [LiveOcean](http://faculty.washington.edu/pmacc/LO/LiveOcean.html).

---

## Why would you clone this repo?

The main users of this repo are people who are in some way collaborating with me and want to use any of my post-processing tools.

---

## Installation

*All the instructions assume you are working from the linux (bash) command line. When I say "go to" I mean navigate to that place, and "do" means enter that command from the linux command line and hit return.*

On your machine, go to wherever you want the LO repo to end up, and do:
```
git clone https://github.com/parkermac/LO.git
```
and LO and all its sub-folders will appear. To get any changes I may make, go to any folder in LO and do:
```
git pull
```
You will also need python3 and various modules. Start by installing the latest version of miniconda on your computer, from [Continuum](https://docs.conda.io/en/latest/miniconda.html). Then create an environment that has all the modules required for this code by going to LO/alpha and executing:
```
conda env create -f loenv.yml
```
You can look in the yml file to see what is being installed.  It even adds the non-python nco toolbox and ffmpeg. For a very few pieces of code you also need matlab, and these are being deprecated, so you can likely ignore this requirement.

---

## Top-level organization

These four directories are assumed to be somewhere, all at the same level in the file structure.

- LO_data: contains large binaries that change infrequently, especially for making grids or forcing files.  I maintain these by hand on my laptop and on my remote linux machines.
- **LO: is this repo.**
- LO_output: is where most output from the LO code ends up, e.g. model forcing files, mooring extractions, plots, etc. It is expected that the contents will change frequently and that they are specific to a given user or machine.
- LO_user: is a placeholder that a user should create to house modified version of the LO code.

LO_output is typically made, if needed, by the code that writes to it. LO_user has to be made by hand (more on that below).

A lot of the code makes use of a dictionary "Ldir" that contains Path objects about where things are. This is created in a somewhat complicated way:
- It is initially specified in `alpha/get_lo_info.py`
- You don't run `get_lo_info.py` itself, but instead it is run every time you run the method `alpha/Lfun.Lstart()` which adds a few more application-specific entries to Ldir.

In principle the code is designed to accommodate other users by allowing them to make another directory (which they can make as their repo) next to LO, called LO_user where you make `alpha/user_get_lo_info.py`.  There are hooks built into the LO code in strategic places that look for
user version of modules allowing you to do things like design your own particle tracking experiments and even add new particle behavior.

#### `get_lo_info.py` or `user_get_lo_info.py` are designed to be the one place where you set machine-dependent choices.  It looks to see what machine you are working on.  It allows you to set several paths to model output, for example: Ldir['roms_out'], Ldir['roms_out1'], and so on.

You can see what is in Ldir by running the module from ipython:
```
run Lfun
```

---

## Organization of LO and relation to LO_output

| LO | LO_output |
| --- | --- |
| alpha: place for things that are shared across folders | |
| pre: pre-processing code, like for loading historical river records | pre/river/[gtag]/... |
| driver: has a couple of drivers that can be used with command-line arguments to (i) create any of the forcing files, and (ii) run one or more ROMS days | |
| forcing: the code for making each of the separate types of forcing | forcing/[gtag]/[fstring]/[frc]/... |
| dot_in: code (one folder for each [gtagex]) for making the .in file for a ROMS run for a given day |
| post: code for automated post-processing of the daily forecast, e.g. for the movies that are sent to the LiveOcean website | post/[gtagex]/layers, etc. |
| extract: code for various types of extractions, plotting, particle tracking, and so on | extract/[gtagex]/cast, etc. |
