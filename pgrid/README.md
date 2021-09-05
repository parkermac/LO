# README for pgrid

### This collection of programs is designed to make gridfiles for ROMS.  It works for analytical and realistic cases, and handles complex tasks like mask editing and smoothing.  It also creates other files associated with rivers, nudging to climatology, and vertical grid parameters, all in the form expected by LO/forcing and ROMS.

---

#### Background on how information about a grid is organized and created.

First edit the gridname and (if needed) a few directory locations at the top of gfun.py.  This gridname will then be used by all subsequent code.  Then most of the programs begin by executing:
```
Gr = gfun.gstart()
```
which returns a dict, Gr, that has the gridname and a few useful Path objects.

In order to keep track of several choices made about a grid, we use dch, a dict of default choices that are initialized using this method:
```
dch =  gfun.default_choices(Gr)
```
You typically override some of the defaults in `start_grid.py` when doing the initial grid specification.  The choices are saved in a pickle file:
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
#### Suggested order to run the code, and bulleted lists of which dch items are used at each step. You have to (1) edit the gridname in `gfun.py`, and (2) edit the grid definition and initial choices in `start_grid.py`.  The rest of the steps are just running the programs in the right order, with `edit_mask.py` requiring user interaction while running.  I would be tempted to make this more automated, but making the grid is a critical part of the model, and many  things can go wrong, so I leave it in individual steps which you can check on along the way using `plot_grid.py`.

`start_grid.py` initializes the grid, bathymetry, and mask. You edit this to define the lon, lat vectors for your grid. It can take a while for big grids, so be patient (e.g. 15 minutes for cas6).
- nudging_edges
- analytical
- t_dir/t_list
- use_z_offset/z_offset

`make_mask.py` makes a first pass at the mask.
- z_land
- unmask_coast
- remove_islands

`edit_mask.py` allows you to fully experience the joy of hand editing a mask.  You run it the first time just to get rid of obvious issues like Lake Washington and river channels.

`carve_rivers.py` uses files of river tracks to make the river channels.

`edit_mask.py` now can be run again for real this time, perhaps running many times. NOTE: you can use an optional command line argument -d ## to change "dval" the carving depth used for lines or points from its default value of 5 m.

`smooth_grid.py` smooths the grid.  This is required for numerical stability.
- use_min_depth/min_depth

`make_extras.py` does some other stuff.
- min_depth (enforced for whole grid)

`grid_to_LO.py` writes the grid.nc file that is used by ROMS, saving it in LO_data/grids. It also saves S-coordinate info.
- nudging_edges
- nudging_days
