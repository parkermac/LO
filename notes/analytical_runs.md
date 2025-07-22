# Analytical (idealized) ROMS simulations

### These are step-by-step directions to guide a new member of the UW Coastal Modeling Group in doing their own ROMS modeling in idealized configurations.

We will be using the LO system of programs, written in python, and you will be working across three separate computers:

1. Your laptop: This is where you will do all your code development, including grid generation.
2. One of our linux machines like perigee or apogee: This is where you will make the forcing files to run ROMS, and where the output will end up.
3. The UW hyak supercomputer klone: This is where the actual ROMS computation is done, making use of the fast parallel nodes we own.

You will be putting your own python installations on some of these, and using GitHub to clone various code repositories to all of them.  It is pretty complicated, but I hope it will be worth the effort. The goal is to arrive at a system that allows you to make effective use of our computing resources, and to make a flexible system for building and analyzing your own ROMS simulations.

### STEPS

#### Install python, the loenv environment, and the LO and LO_user repos

Install python and the LO code repo on your laptop. Even if you already have a python installation this should not interfere with what you have because we will use an environment. Follow the directions [HERE](https://github.com/parkermac/LO/blob/main/README.md). Now you have the "loenv" python environment that has access to all the python modules you need, and you have the LO and LO_user code repos.

Note that you will have to create your own repo called `LO_user`, then copy `get_lo_info.py` into it, and edit it with information specific to you. This is the one place where user-specific and machine-specific information is kept.

#### Make the ROMS grid for an analytical run

Now make a grid using the pgrid programs in LO. First make a directory in your LO_user called pgrid, and then copy `LO/pgrid/gfun_user.py` into `LO_user/pgrid/gfun_user.py` and check that the "ae0" grid is the one it is working on. Around line 19 you should see:
```
gridname = 'ae0'
```
and if not, then edit it to make it so. This is an idealized estuary + coast domain that looks like this: ![ae0](./figures/ae0.png) You can look farther down in the code and see how this is generated.  Note that we also create a little pandas DataFrame with info about a river called "creek0", and some track info.

To make this grid yourself, go to `LO/pgrid`, launch ipython, and run these commands in order:
- start_grid
- make_mask
- carve_rivers
- smooth_grid
- make_extras
- grid_to_LO

And using plot_grid you should be able to generate the figure above.  When you run each command you are presented with a list of possible grid iterations to work on; by hitting return each time you choose the last iteration, which is almost always what you want to do. There is a detailed README in pgrid to help with what each step is doing. One thing we did not do in the case was run edit_mask.  For realistic runs getting the mask to match the coastline needs some user attention, but in this case it is not needed.

You may notice along the way that the code made **new directories** on your computer: `LO_data/grids/ae0` and `LO_output/pgrid`. These are instances of the general [organizing structure](http://faculty.washington.edu/pmacc/Research/new_ideas.html) I use for all projects:

[Data] -> [Code] -> [Output]

---

**NOTE: the handling of river information can be confusing. Here is what happens:**

1. LO/pgrid/start_grid.py (calling gfun_user.py) creates files that have information about rivers (names and gage numbers - if any) and their channel locations. For non-analytical runs these are would be created by the LO/pre/river1 programs. The README in LO/pre/river1 gives a lot more information. For analytical cases things are much simpler so we just make these by hand using code in gfun_user.py.
    - LO_output/pre/river1/ae0_v0/river_info.csv
    - LO_output/pre/river1/ae0_v0/tracks/creek0.p
2. LO/pgrid/carve_rivers.py goes through each river and figures out the exact grid indices where its source should be on the ROMS grid you are making, saving the output in:
    - LO_output/pgrid/ae0/roms_river_info.csv
3. LO/pgrid/grid_to_LO.py copies the roms_river_info.csv you just made into a place where it can be used by the forcing code, e.g. rivA0.
    - LO_data/grids/ae0/river_info.csv

The confusing bit is that the name "river_info.csv" was reused in steps 1 and 3, even though in 3 it is a file full of indices, and in 1 it has gage numbers. This is just a legacy of how the code was developed, and me not wanting to break things already in place.

---

#### Repeat some steps on a remote linux machine

Get an account on perigee or apogee from David Darr.  Then, working in the main directory he makes for you (e.g. /dat1/[username], not in your $HOME directory ~) again install python and the LO (clone from parkermac) and LO_user (clone from your GitHub account) repos.

Then copy LO_data/grids/ae0 from your laptop to the remote machine.  You could use scp for this, or I like the Transmit program from the App store, but it is $25 per year.

---

#### Create the forcing files for a run

Working on a remote machine like apogee or perigee, go to `LO/driver` and execute these three commands from the linux command line (sequence doesn't matter):
```
python driver_forcing3.py -g ae0 -0 2020.01.01 -1 2020.01.02 -f rivA0
python driver_forcing3.py -g ae0 -0 2020.01.01 -1 2020.01.02 -f ocnA0 -s new
python driver_forcing3.py -g ae0 -0 2020.01.01 -1 2020.01.02 -f tideA0
```
These will make files in `LO_output/forcing/ae0/f2020.01.01` and `f2020.01.02`. These are the NetCDF files ROMS will need to force the run for two days.

You can also run these on your laptop. The reason you need to run them on the remote machine is that klone will be looking on the remote machine when it tries to get forcing for a given day.

#### Run ROMS on klone

See [HERE](https://github.com/parkermac/LO_roms_user/blob/main/README.md) for detailed info on getting ROMS working on klone. You should follow these all the way through running the upwelling example and creating your own LO_roms_user repo.

Compile the xa0 executable.

Copy the grid file to klone with a command like this:
```
scp -r [username]@apogee.ocean.washington.edu:/dat1/[username]/LO_data/grids/ae0 .
```
Run the analytical case, from the head node, with a command like this.
```
python3 driver_roms4.py -g ae0 -t t0 -x xa0 -s newcontinuation -0 2020.01.01 -1 2020.01.02 --group_choice macc --cpu_choice compute -np 40 -N 40 < /dev/null > ae.log &
```
