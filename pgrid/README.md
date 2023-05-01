# README for pgrid

### This collection of programs is designed to make gridfiles for ROMS.  It works for analytical and realistic cases, and handles complex tasks like mask editing and smoothing.  It also creates other files associated with rivers, nudging to climatology, and vertical grid parameters, all in the form expected by LO/forcing and ROMS.

Throughout the steps we maintain and modify the arrays "h" and "mask_rho" to be complete 2-D arrays with no nan's.  mask_rho follows the ROMS convention of having 1 = water, 0 = land.

To see what variables ROMS expects to find in a grid file, have a look in `LO_roms_source/ROMS/Utility/get_grid.F`.

NOTE: As of 2023.04.13 I have moved over to the new river data format in LO/pre/river1. This meant introducing a new dch['ctag'] which corresponds to the collection tag in river1. The current default is lo_base. We also replace Gr['ri_dir'] with Gr['ri_dir0'] which is the path to LO_output/pre/river1.

---
#### Suggested order to run the code

NOTE: There is a streamlined set of these instructions in `LO/notes/analytical_runs.md` that are specific to getting an analytical run going.  You might want to follow them the first time, and then refer to these instructions when you get deeper into making your own grids.

(0) Start by making your own LO_user repo (see `LO/README.md`) and then adding a sub-folder "pgrid", and then copying `LO/pgrid/gfun_user.py` into it.

(1) You have to edit the gridname and definition in `LO_user/gfun_user.py`.

There are a number of choices you can set as parameters.  You can see these in the module `gfun.default_choices()`.  They end up in a dict called dch ("default choices") that gets saved and reopened at each step.  You can override any of the defaults in your `gfun_user.py` entry.

Note: you make your grid definition as an if-statement in `gfun_user.py` AND in the line near the top of that module "gridname = ".

The rest of the steps are just running the programs in the right order, with `edit_mask.py` requiring user interaction while running.  I would be tempted to make this more automated, but making the grid is a critical part of the model, and many  things can go wrong, so I leave it in individual steps which you can check on along the way using `plot_grid.py`.

Throughout this code I try to use ROMS naming conventions, except that when manipulating or plotting I refer to: [lon_rho, lat_rho] as [lon, lat].

Also [plon, plat] is just like [lon_psi, lat_psi] but extended by one on all directions so that it is box corners around all rho-grid points.

- The bulleted lists below each step are the dch items used in that step.

(2) `start_grid.py` initializes the grid, bathymetry, and mask. You edit this to define the lon, lat vectors for your grid.
- nudging_edges
- analytical
- t_dir/t_list
- use_z_offset/z_offset
- excluded_rivers
- any other of the default choices in `gfun.default_choices()`

(3) `make_mask.py` makes a first pass at the mask.
- z_land
- unmask_coast
- remove_islands

(4) `edit_mask.py` allows you to fully experience the joy of hand editing a mask.  You run it the first time just to get rid of obvious issues like Lake Washington and river channels.

(5) `carve_rivers.py` uses files of river tracks to make the river channels.

You can set a list of rivers to exclude in `gfun_user.py`.

There is a new tool called `create_river_tracks.py` which allows you to define, or redefine the river tracks.

(6) `edit_mask.py` now can be run again for real this time, perhaps running many times. NOTE: you can use an optional command line argument -d ## to change "dval" the carving depth used for lines or points from its default value of 5 m.

(7) `smooth_grid.py` smooths the grid.  This is required for numerical stability.
- use_min_depth/min_depth

NOTE: you may want to run `carve_rivers.py` again at this point just to make just they are still there.

(8) `make_extras.py` makes the u, v, and psi masks.  Note: you cannot do any mask editing after this step, because it creates the other masks based on mask_rho.  If you do edit the mask, run this again.
- use_min_depth/min_depth (enforced for whole grid)

(9) `grid_to_LO.py` writes the grid.nc file that is used by ROMS, saving it in LO_data/grids. It also saves S-coordinate info.

It copies or creates: grid.nc, nudgcoef.nc, river_info.csv, S_COORDINATE_INFO.csv, XY_COORDINATE_INFO.csv. NOTE: the river_info.csv file contains the index and direction info used by ROMS. Don't confuse it with the file of the same name (sorry) that is in LO_output/pre/river.
- nudging_edges
- nudging_days

---

#### Background on how information about a grid is organized and created.

Most of the programs begin by executing:
```
Gr = gfun.gstart()
```
which returns a dict, Gr, that has the gridname and a few useful Path objects.

In order to keep track of several choices made about a grid, we use dch, a dict of default choices that are initialized using this method:
```
dch =  gfun.default_choices()
```
You typically override some of the defaults in `LO_user/pgrid/gfun_user.py` when doing the initial grid specification.  The choices are saved in a pickle file:
```
pickle.dump(dch, open(Gr['gdir'] + 'choices.p', 'wb'))
```
You can also go back and change things in dch later using `Z_edit_dch.py`.

Each time you run a piece of code it makes a new grid file with the name altered to indicate what happened.  The names start as: grid_m00_r00_s00_x00.nc with the letters and numbers indicating changes to: mask (m), river (r), smoothing (s), or extras (x).

You can use `plot_grid.py` to look at any of the grids.  You can override the default grid by calling it with a flag:
```
run plot_grid.py -g sal0
```

---

#### Development Notes

Throughout this code I use xr.open_dataset() and then the .update method to store changed variables before saving to a new file.  This works because I am saving to a new NetCDF file anytime there are changes. In xarray if you are changing data in a file and saving to the same name you have to use the more memory-intensive xr.load_dataset().
