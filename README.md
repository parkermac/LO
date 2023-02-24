# README for LO

## LO is the code base for running the LiveOcean collection of regional ocean simulations.

The code here handles all the non-research tasks of making model (e.g. ROMS) forcing files, running the model on a remote linux cluster, and post-processing of model output.  Many of the post-processing tasks, such as mooring- cast- and layer extractions, and particle tracking, are designed for other users to be able to run from the linux command line, or modify if needed.

This version is the second major release of this project, the first was called LiveOcean. This version incorporates all the things I learned building the first one. A few specifics: it uses python Path objects everywhere possible. It defaults to using NetCDF4 for all the ROMS forcing. It uses python code instead of shell scripts wherever that makes things simpler - especially in the drivers. It also uses a much more rigorous and naming and organization system. It transitions to using xarray instead of the netCDF4 module.  For a very few pieces of code you also need Matlab, and these are being deprecated, so you can likely ignore this requirement.

All the instructions assume you are working from the linux (bash) command line.

Contact Parker MacCready, p.maccready@gmail.com, with any questions.

Check out today's model output at [LiveOcean](http://faculty.washington.edu/pmacc/LO/LiveOcean.html).

---

## Why would you clone this repo?

The main users of this repo are people who are in some way collaborating with me and want to use any of my forcing-generation, model-running, or post-processing tools. You don't need this repo to work with the LiveOcean ROMS output files, but some of the tools here may be useful.

---

## Getting up to speed with linux, python, and git

To work with the LO system, ROMS, and generally to be able to do modern numerical modeling, it is important that you have reasonable skill with linux, python, and git.  If you need a refresher or a starting place with any of these you could look at the lectures from my [Effective Computing Class](http://faculty.washington.edu/pmacc/Classes/EffCom_2020/index.html).  Ignore the first lectures about installing linux and python.  They are superseded by the notes here.  But the lectures starting with [Linux 1](http://faculty.washington.edu/pmacc/Classes/EffCom_2020/lectures/Linux%201.pdf) and continuing through those on python and GitHub will be relevant.

---

## Installation - four steps

All the instructions assume you are working from the linux (bash) command line. When I say "go to" I mean navigate to that place, and "do" means enter that command from the linux command line and hit return.

For mac users you already have the linux operating system and a terminal. To help with cross-platform compatibility we will used the bash shell (instead of the new mac default zsh), so you should issue this command in your terminal window `chsh -s /bin/bash` to set your default shell to bash.

Windows users need to download Ubuntu and get it configured. Here is a [LINK](https://ubuntu.com/download/desktop) to the download. After downloading, follow these [DIRECTIONS](https://docs.microsoft.com/en-us/windows/wsl/install) about the best practice to get this set up before installing miniconda. WSL is basically configuring linux within windows for ubuntu to use.  Then just work from a terminal window using ubuntu as your version of linux. _Thanks to Marissa Leatherman for these tips._

To use the LO system effectively as a member of MacCready's research group, you will likely be setting up your own python installation on both your laptop and on one or more of our group servers (e.g. apogee and perigee).  You will also probably be working on the UW supercomputer cluster called hyak, but you don't need to install python there.

#### (1) Install miniconda on your laptop and then create the "loenv" environment.

Anaconda is a great way to get python. We will use a minimal installation of python3 called "miniconda" and then create a conda "environment" and add required packages to it ourselves.  This keeps things lightweight and less likely to suffer from package conflicts.

Go to this miniconda webpage to find the link for the "latest" installer for your platform: [INSTALLERS](https://docs.conda.io/en/latest/miniconda.html#latest-miniconda-installer-links).  The "latest" is the table at the top of the page.  The bash versions are shell scripts, and the pkg versions are graphical installers.  These directions assume you are using the bash shellscript.

On a Mac click on the link to download the installer, or copy the link address and then from the command line do:
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
and make sure that the number matches the one from the table of installers.  Note: this may only work on linux machines; I don't have sha256sum on my mac.

Next run the installer:
```
bash [installer.sh]
```
and answer yes to everything, especially "initialization".  For one of my linux installations I chose to put it in /data1/parker/miniconda3, or on my mac I put it in /Applications/miniconda3.  It doesn't matter what folder you put it in.  The initialization adds some lines to your ~/.bashrc or ~/.bash_profile (mac). After you "source" this you should have (base) appended to your bash prompt, and this tells you that you are in the (base) conda environment so you can run python - but not much else, yet.

Then for good measure do:
```
conda update conda
```

You will also need git, which you may or may not have.  To find out if you have it, type `which git` in the terminal and see if it finds anything.  If you need to get it you can install it easily by doing:
```
conda install -c anaconda git
```

#### (2) Get the LO code

The main pile of programs that you will use are in a GitHub repo called LO, maintained by MacCready.  You will clone this to all the machines you are working on, and should not have to make any changes to it, except for updating occasionally.

On your machine, first go to wherever you want the LO repo to end up.  On my mac I go to Documents.  On one of our linux machines you might put it in /data1/[username] which is the working directory we made for you.  On the linux mashines don't put it in your home (~) folder because there is not much disk space there. When you are at whatever place you have decided on do:
```
git clone https://github.com/parkermac/LO.git
```
and LO and all its sub-folders will appear. To get any changes I may make, go to any folder in LO and do:
```
git pull
```
Sometimes the remote machines will complain about this. One solution is to issue this command while inside the repo on the remote machine:
```
git config pull.ff only
```
which (I think) means that your pull always does "fast forward" moving the local branch to the same state as the one in the cloud. This is appropriate if you are using the "one way" git workflow that I generally follow: only edit code on my laptop, push to the cloud, pull from the cloud to remote machines; **never edit on remote machines**.

#### (3) Create the (loenv) environment

Create an environment that has all the modules required for this code by going to LO and executing:
```
conda env create -f loenv.yml > env.log &
```
This may take a half hour or so, which is why I added the `> env.log &` at the end of the command. You can look in the .yml file to see what is being installed.  It even adds the non-python nco toolbox and ffmpeg.  It also adds LO/lo_tools as a local "package" so that when you are in (loenv) you can access any of the modules in LO/lo_tools/lo_tools with a line in your python code like `from lo_tools import zfun`.  Instructions for a simple approach to making your own local packages can be found [HERE](https://pythonchb.github.io/PythonTopics/where_to_put_your_code.html).

If you like you can make your own .yml and make your own environment, especially if you want to add additional packages.  The LO code does not need to be run in (loenv) but it does assume that lo_tools is an installed local package. One question I have is, if you made your own `LO_user/loenv.yml` file you need to set a path in the pip line that installs the lo_tools package (see instructions below, under: Creating your own environment).

Then if you want to use this environment all the time add this line:
```
conda activate loenv
```
to your .bashrc or .bash_profile, and "source" it.  Now (loenv) will appear at the start of your bash prompt.

#### _Updating the environment_

Generally whenever I am about to update a python environment, I first update conda, doing:
```
conda update -n base -c defaults conda
```
Then to update the loenv environment, go to LO where the loenv.yml file is and do:
```
conda env update -f loenv.yml
```
which only took a minute or two the last time I tried it. There are more complete instructions for working with conda environments [HERE](https://www.earthdatascience.org/courses/intro-to-earth-data-science/python-code-fundamentals/use-python-packages/use-conda-environments-and-install-packages/).

NOTE: on 2023.01.13 I encountered a problem where pickled pandas DataFrames created with pandas 1.2 and earlier could not be opened with pandas 1.3 and higher. For some reason when I updated everything on perigee the pandas version remained at 1.2, even though it was 1.3 on apogee and my mac. To solve this I did a hand install on perigee:
```
pip install --upgrade pandas=1.3.3
```
and this fixed the problem without breaking anything that I noticed. Hopefully this is a temporary inconsistency.

NOTE: 2023.02.24 Kate suggests that using pip3 in place of pip (including in the .yml file) may be more up-to-date (always uses python 3 version). I have not tested this yet.

#### _Creating your own environment! Highly recommended_

One way to do this would be to:

- copy `loenv.yml` to LO_user or any other repo that is yours
- rename it, for example, to `myenv.yml`
- edit it so that `name: myenv` is the first line
- add or subtract any packages you like
- if you are using the lo_tools local package, change the path for it in the yml to be `-e ../LO/lo_tools lo_tools` (assuming that is the correct relative path)
- then do: `conda env create -f myenv.yml > env.log &` and there you are!

At any time you can do `conda info --envs` to find out what environments you have. And if you want to cleanly get rid of an environment, just make sure it is not active, then do `conda env remove -n myenv` or whatever name you are wanting to remove. This will not delete your yml file.

#### _NOTE for mac users_

As of 2022.12.07 the conda-forge version of the nco operators does not have a version for the new M1 chip. So you will want to make your own environment yaml file, deleting nco from the list (and pytide as well, with is not used much). And then install nco using homebrew instead:
- Instructions to install homebrew: https://brew.sh/
- Install nco using homebrew: https://formulae.brew.sh/formula/nco#default
- Both are one-liners. Easy!

After installing nco using homebrew, you can add -nco back to your myenv.yml and in terminal do: conda env update -f myenv.yml. A similar solution might exist for pytide, but I haven't tried it yet.

#### (4) Create your own LO_user and make it a GitHub repo

For this you will be working on your laptop and creating your own GitHub repo.  I can't stress enough how useful it is to keep all your code in git. It saves your code to the cloud. It allows you to see changes you have made. It makes it easy to share code with others.  Finally it makes it easy to maintain your code on other machines.

Here are some [INSTRUCTIONS](http://faculty.washington.edu/pmacc/Classes/EffCom_2020/lectures/GitHub%20Intro.pdf) for how to get started using Git and making your own repo.  These are for using the GitHub desktop app.  You can also accomplish all the same things using command line statements following these [INSTRUCTIONS](https://docs.github.com/en/get-started/importing-your-projects-to-github/importing-source-code-to-github/adding-locally-hosted-code-to-github).

**In either case, you will need to create your own account on GitHub.**

LO_user is a place where the LO code looks for user versions of things, like particle tracking experiment initial conditions.  Most importantly, it is where user- and machine-specific paths are defined in `get_lo_info.py`.

**First**: at the same level as LO on your laptop, create the directory LO_user.

**Second**: using either GitHub Desktop or the command line, make LO_user into a git repo on you laptop (initialize and commit).  I usually start a repo in GitHub Desktop by choosing the "python" .gitignore from the defaults, along with the MIT license, and making the repo public.  Then publish to your account on GitHub in the cloud.

PRO TIP: If you are working on a mac you will find that the .DS_Store file clutters up your repo (it is just saving folder view choice information).  To de-clutter this:
- Add .DS_Store to your .gitignore.
- If you have instances of .DS_Store already in your repo, you can clean them out from the command line by going to the repo and doing:
```
find . -name .DS_Store -print0 | xargs -0 git rm --ignore-unmatch -f
```

**Third**: Copy `get_lo_info.py` from LO to LO_user, and make a couple of edits to these lines:
```
if str(HOME) == '/Users/pm8':
    lo_env = 'pm_mac'
```
- change `/Users/pm8` to whatever you get by doing `echo $HOME` from your linux command line.
- change `pm_mac` to some string that is [your initials or whatever]_[mac or pc].  The mac or pc part is important because there is a place in the LO plotting code where it looks for these and makes a decision about the graphical backend.

Check that things are working as you expect by going to LO/driver and doing:
```
python test_loenv.py
```
You should get some screen output of the contents of the python dict called Ldir.  If this does not work, email Parker for help.

If everything is working, go ahead and push this to your git repo in the cloud.

---

#### Hooks connecting LO and LO_user

The user version of `get_lo_info.py` that you created above in LO_user is an example of a "hook" built into the LO code.  In this case the `lo_tools/Lfun.Lstart()` method that is called at the top of most LO code looks first to see if `LO_user/get_lo_info.py` exists, and if so it uses it to fill to fill out the Ldir dict.  Otherwise it uses the default (Parker's) `LO/get_lo_info.py`.

Similar hooks are built into other parts of the LO system where we expect that users will want to use the LO code but have it access customized versions of some parts.  Here is the list of current hooks that it looks for:
- `LO_user/get_lo_info.py`
- `LO_user/pgrid/gfun_user.py`
- `LO_user/forcing/[frc]` (used by `LO/driver_forcing.py`)
- `LO_user/dot_in/[gtagex]` (used by `LO/driver_roms2.py`)
- `LO_user/tracker/experiments.py` and `LO_user/tracker/trackfun.py`
- `LO_user/extract/box/job_definitions.py` (used by `LO/extract/box/extract_box.py` and `extract_box_chunks.py`)
- `LO_user/extract/moor/job_lists.py` (used by `LO/extract/moor/multi_mooring_driver.py`)
- `LO_user/plotting/roms_plots.py` (used by `LO/plotting/pan_plot.py`)

Not yet implemented, but should be:
- Ask for what you want to see!

NOTE: there is now a [blog-like file](https://github.com/parkermac/LO/blob/main/notes/USER_NOTES.md) `LO/notes/USER_NOTES.md` where I keep notes about changes to the code that users should be aware of.

NOTE: The hooks are meant to allow a user to change a single part of the LO system (like having their own lists of mooring locations) while still making use of the generic and well-tested LO machinery (like for mooring extraction). But for other types of analysis the user may be better off just copying a whole directory of code into LO and then editing it to do what they need. For example, the `LO/extract/tef` code is probably too closely tied to the TEF sections and segments from the MacCready et al. (2021, JGR) paper. If you do make a copy of an LO code folder and put it in LO_user, I suggest that you maintain the LO directory structure.  In the example this would mean putting the code in `LO_user/extract/tef`. Most LO code is designed to load modules either from lo_tools (which you always have access to if you are in the loenv python environment) or from modules in the current working directory. For output, unless you change things it will probably go to LO_output, just like the things from LO. I think this is probably fine - no need to create a separate LO_user_output.

---

## Top-level organization

These four directories are assumed to be somewhere, all at the same level in the file structure.

- LO: is this repo.
- LO_user: is a required separate folder for information and programs specific to a given user.
- LO_data: contains large binaries that change infrequently, especially for making grids or forcing files.  I maintain these by hand on my laptop and on my remote linux machines. It may be advantageous for a user on those remote machines to just point to my LO_data (see `get_lo_info.py`) instead of copying to make their own.  Of course on their personal laptop they will need to make their own LO_data.
- LO_output: is where most output from the LO code ends up, e.g. model forcing files, mooring extractions, plots, etc. It is expected that the contents will change frequently and that they are specific to a given user or machine. LO_output is typically made, if needed, by the code that writes to it.

A lot of the code makes use of a dictionary "Ldir" that contains user- and machine-specific Path objects about where things are. This is created in a somewhat complicated way:
- It is initially specified in `LO_user/get_lo_info.py`
- You don't run `get_lo_info.py` itself, but instead it is run every time you run the method `lo_tools/lo_tools/Lfun.Lstart()` which adds a few more application-specific entries to Ldir.

`get_lo_info.py` is designed to be the **one place** where you set machine-dependent choices.  It looks to see what machine you are working on.  It allows you to set several paths to model output, for example: Ldir['roms_out'], Ldir['roms_out1'], and so on.

---

## Guide to the various READMEs

Most of the sub-folders of LO have their own README files, but a few are useful to be aware of from the very start:

- `LO/README.md` (this file) Initial python installation, cloning of the LO code, and creation of your LO_user repo.
- `LO/notes/analytical_runs.md` Step-by-step instructions for making and running your own analytical (idealized) ROMS run. This covers grid generation on your laptop using LO/pgrid. and the creation of forcing files on one of our servers using various bits in LO/forcing.  Then it covers installing ROMS on hyak, testing that it works, and compiling it and running it for your analytical grid.
- `LO/notes/ROMS_Tips.md` has some useful hints about finding your way around the important parts of the ROMS source code.
- `LO/notes/Forecast_Operators_Manual.md` has instructions for checking on and troubleshooting the LiveOcean daily forecast.
- `LO_roms_user/README.md` Notes about setting up your environment on hyak, installing ROMS, and running the upwelling test case.  It also gives a listing of LO compiler configurations in Parker's repo.

---

## What to put where?

Running the LO code on our model output or running ROMS yourself require that you work across two or three computers (laptop : perigee/apogee : mox/klone) it is useful to know what things you need to put on which machines.  Here is an outline:

First, do the miniconda/loenv python installation on both your laptop and on perigee/apogee. You do not need to do this on mox/klone because python is already installed there and we don't require extra packages.

In the fields below:
- 1 = your laptop
- 2 = perigee or apogee
- 3 = mox or klone

| Folder | How created | Put where | Purpose |
| --- | --- | --- | --- |
| LO   | Clone from parkermac | 1,2,3  | Primary tools for LiveOcean  |
| loenv   | See instructions above - each user creates their own from scratch  | 1,2  | python environment  |
| LO_user   | Create your own repo on 1 and clone to 2,3, then copy in bits from LO and edit them.  | 1,2,3  | User versions of LO code  |
| LO_data   | Create your own folder on 1 and copy to 2,3  | 1,2,3  | Data files like grids and coastlines. Not a repo.  |
| LO_output   | Created automatically  | Will appear on 1,2,3  | LO output, like forcing files. Not a repo.  |
| LO_roms   | Created automatically  | Will appear on 2,3 and you can copy to 1  | ROMS output. Not a repo.  |
| LO_roms_source   | Clone using svn from ROMS site.   | 1,3  | ROMS source code. Good to have on 1 for reference, even though you only compile on 3  |
| LO_roms_source_alt   | Clone from parkermac  | 1,2,3  | Our customized bits of the ROMS source code, especially the biology routines.  Also an edited version of varinfo.yaml  |
| LO_roms_user   | Create your own repo on 1 and copy in bits from the verion on parkermac. Then clone to 3.  | 1,3  | Files for configuring specific ROMS versions that are defined by lists of compiler flags.  |

---

## Organization of LO and relation to LO_output

**Naming convention**: We follow a strict system of naming things associated with a ROMS run in order to allow for modularity, e.g. running a given grid and forcing files with a new executable.

Things that I type in [ ] below mean that they would be replaced by specific strings, for example when using them as command line arguments.

- [gridname] is the name of the grid (e.g. cas6)
- [tag] is a name to identify any of the things that are controlled by a "dot_in" instance (e.g. v00)
- [ex_name] is the name of the ROMS executable (e.g. u0kb)
- [fstring] is a date string of the form fYYYY.MM.DD (e.g. f2021.07.04)
- [date_string] is a date string of the form YYYY.MM.DD (e.g. 2021.07.04)
- [frc] is the name of one of the model forcings (e.g. ocn0)

Grids are just identified by [gridname].

Collections of forcing files are identified by [gridname] since they are always created for a specific grid. Note that you can accumulate many different types of forcing inside a certain [gridname], e.g. ocn00, ocn01, etc. This naming logic was introduced with `driver_forcing3.py` and the associated `driver_roms3.py`.

A specific run is identified by [gridname]\_[tag]\_[ex_name] which is also referred to as [gtagex] or Ldir['gtagex']. The [tag] is only introduced in the creation of a dot_in instance, which is in a folder named LO/dot_in/[gtagex]. The [tag] represents the specific set of forcing files and other dot_in choices that define a run.

NOTE: some of the code will parse a gtagex into it constituent gridname, tag, and ex_name.  To do so it assumes these are separated by an underscore "_", so don't use underscores in any of your gridnames, tags, or ex_names.

Here is some info on the various folders in LO, and how they relate to the naming of things that end up in LO_output.  Many of these have their own README files.

| LO | LO_output |
| --- | --- |
| lo_tools/lo_tools: place for shared modules | |
| pre: pre-processing code, like for loading historical river records | pre/river/cas6_v3/... |
| driver: has a couple of drivers that can be used with command-line arguments to (i) create any of the forcing files, and (ii) run one or more ROMS days | |
| forcing: the code for making each of the separate types of forcing | forcing/[gridname]/[fstring]/[frc]/... |
| dot_in: code (one folder for each [gtagex]) for making the .in file for a ROMS run for a given day |
| post: code for automated post-processing of the daily forecast, e.g. for the movies that are sent to the LiveOcean website | post/[gtagex]/layers, etc. |
| extract: code for various types of extractions | extract/[gtagex]/cast, etc. |
| pgrid: code for making ROMS grids   |   |
| tracker: particle tracking code   |   |
| daymovie: code to make the daily forecast movies for the website   |   |
| notes: README's on various topics   |   |
