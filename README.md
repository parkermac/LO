# LO is the code base for running the LiveOcean collection of regional ocean simulations.

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


---

## Top-level organization

These four directories are assumed to be somewhere, all at the same level in the file structure.

- LO_data: contains large binaries that change infrequently, especially for making grids or forcing files.  I maintain these by hand on my laptop and on my remote linux machines.
- **LO: is this repo.**
- LO_output: is where most output from the LO code ends up, e.g. model forcing files, mooring extractions, plots, etc. It is expected that the contents will change frequently and that they are specific to a given user or machine.
- LO_user: is a placeholder that a user should create to house modified version of the LO code.

LO_output is typically made, if needed, by the code that writes to it. LO_user has to be made by hand (more on that below).

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
