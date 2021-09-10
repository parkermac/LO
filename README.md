# README for LO

## LO is the code base for running the LiveOcean collection of regional ocean simulations.

The code here handles all the non-research tasks of making model (e.g. ROMS) forcing files, running the model on a remote linux cluster, and post-processing of model output.  Many of the post-processing tasks, such as mooring- cast- and layer extractions, and particle tracking, are designed for other users to be able to run from the linux command line, or modify if needed.

This version is the second major release of this project, the first was called LiveOcean. This version incorporates all the things I learned building the first one. A few specifics: it uses python Path objects everywhere possible. It defaults to using NetCDF4 for all the ROMS forcing. It uses python code instead of shell scripts wherever that makes things simpler - especially in the drivers. It also uses a much more rigorous and naming and organization system. It transitions to using xarray instead of the netCDF4 module.  For a very few pieces of code you also need Matlab, and these are being deprecated, so you can likely ignore this requirement.

All the instructions assume you are working from the linux (bash) command line.

Contact Parker MacCready, p.maccready@gmail.com, with any questions.

Check out today's model output at [LiveOcean](http://faculty.washington.edu/pmacc/LO/LiveOcean.html).

---

## Why would you clone this repo?

The main users of this repo are people who are in some way collaborating with me and want to use any of my forcing-generation, model-running, or post-processing tools.

---

## Installation - four steps

*All the instructions assume you are working from the linux (bash) command line. When I say "go to" I mean navigate to that place, and "do" means enter that command from the linux command line and hit return.*

#### First: install miniconda on your machine and then create the "loenv" environment.

Anaconda is a great way to get python. We will use a minimal installation of python3 called "miniconda" and then create a conda "environment" and add required packages to it ourselves.  This keeps things lightweight and less likely to suffer from package conflicts.

Go to this miniconda webpage to find the link for the "latest" installer for your platform: [INSTALLERS](https://docs.conda.io/en/latest/miniconda.html#latest-miniconda-installer-links).  The "latest" is the table at the top of the page.  For example, for linux it is currently: https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh.

On a Mac click on the link to download the installer or copy the link address and from the command line do:
```
curl -O [address]
```
Or, from the linux command line do:
```
wget [address]
```
This will get you a file that I'll call [installer.sh].  It doesn't really matter where you put the installer.  You can choose where the installation will end up when you run it.

Then if you are being really careful about security do:
```
sha256sum [installer.sh]
```
and make sure that the number matches the one from the table of installers.

Next run the installer:
```
bash [installer.sh]
```
and answer yes to everything, especially "initialization".  For one of my linux installations I chose to put it in /data1/parker/miniconda3, or on my Mac I put it in /Applications/miniconda3.  The initialization adds some lines to your ~/.bashrc or ~/.bash_profile (Mac). After you "source" this you should have (base) appended to your bash prompt, and this tells you that you are in the (base) conda environment so you can run python - but not much else, yet.

Then for good measure do:
```
conda update conda
```

You will also need git, which you may or may not have.  You can install it easily by doing:
```
conda install -c anaconda git
```

#### Second: Get the LO code

On your machine, go to wherever you want the LO repo to end up, and do:
```
git clone https://github.com/parkermac/LO.git
```
and LO and all its sub-folders will appear. To get any changes I may make, go to any folder in LO and do:
```
git pull
```

#### Third: create the (loenv) environment

Create an environment that has all the modules required for this code by going to LO and executing:
```
conda env create -f loenv.yml > env.log &
```
This may take a half hour or so, which is why I added the `> env.log &` at the end of the command (I haven't actually tested that yet). You can look in the .yml file to see what is being installed.  It even adds the non-python nco toolbox and ffmpeg.  It also adds LO/lo_tools as a local "package" so that when you are in (loenv) you can access any of the modules in LO/lo_tools/lo_tools with a line in your python code like `from lo_tools import zfun`.  Instructions for a simple approach to making your own local packages can be found [HERE](https://pythonchb.github.io/PythonTopics/where_to_put_your_code.html).

If you like you can make your own .yml and make your own environment, especially if you wan to add additional packages.  The LO code does not need to be run in (loenv) but it does assume that lo_tools is an installed local package.

Then if you want to use this environment all the time add this line:
```
conda activate loenv
```
to your .bashrc or .bash_profile, and "source" it.  Now (loenv) will appear at the start of your bash prompt.

#### Fourth: create LO_user

Make a folder called LO_user at the same level as LO, and then put a file inside it called get_lo_info.py. You can copy this from my own LO_user repo: [FILE](https://github.com/parkermac/LO_user/blob/main/get_lo_info.py).

Here are some [INSTRUCTIONS](http://faculty.washington.edu/pmacc/Classes/EffCom_2020/lectures/GitHub%20Intro.pdf) for how to get started using Git and making your own repo. **I strongly encourage you to make LO_user your own repo in GitHub.**

LO_user will also be a place where the LO code looks for user versions of things, like particle tracking experiment initial conditions.

---

## Top-level organization

These four directories are assumed to be somewhere, all at the same level in the file structure.

- LO_data: contains large binaries that change infrequently, especially for making grids or forcing files.  I maintain these by hand on my laptop and on my remote linux machines.
- **LO: is this repo.**
- LO_output: is where most output from the LO code ends up, e.g. model forcing files, mooring extractions, plots, etc. It is expected that the contents will change frequently and that they are specific to a given user or machine.
- LO_user: is a required separate folder for information and programs specific to a given user.

LO_output is typically made, if needed, by the code that writes to it. LO_user has to be made by hand (more on that below).

A lot of the code makes use of a dictionary "Ldir" that contains Path objects about where things are. This is created in a somewhat complicated way:
- It is initially specified in `LO_user/get_lo_info.py`
- You don't run `get_lo_info.py` itself, but instead it is run every time you run the method `lo_tools/lo_tools/Lfun.Lstart()` which adds a few more application-specific entries to Ldir.

#### `get_lo_info.py` is designed to be the one place where you set machine-dependent choices.  It looks to see what machine you are working on.  It allows you to set several paths to model output, for example: Ldir['roms_out'], Ldir['roms_out1'], and so on.

---

## Organization of LO and relation to LO_output

To Do: need to explain what [gtag] etc. are.

| LO | LO_output |
| --- | --- |
| lo_tools/lo_tools: place for shared modules | |
| pre: pre-processing code, like for loading historical river records | pre/river/[gtag]/... |
| driver: has a couple of drivers that can be used with command-line arguments to (i) create any of the forcing files, and (ii) run one or more ROMS days | |
| forcing: the code for making each of the separate types of forcing | forcing/[gtag]/[fstring]/[frc]/... |
| dot_in: code (one folder for each [gtagex]) for making the .in file for a ROMS run for a given day |
| post: code for automated post-processing of the daily forecast, e.g. for the movies that are sent to the LiveOcean website | post/[gtagex]/layers, etc. |
| extract: code for various types of extractions, plotting, particle tracking, and so on | extract/[gtagex]/cast, etc. |
