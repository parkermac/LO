This collection of programs is designed to make gridfiles for ROMS.  It works for analytical and realistic cases, and handles complex tasks like mask editing and smoothing.  It also creates other files associated with rivers, nudging to climatology, and vertical grid parameters, all in the form expected by LiveOcean/forcing and ROMS.

First edit the gridname and (if needed) a few directory locations at the top of gfun.py.  This gridname will then be used by all subsequent code.

In order to keep track of several choices typically made about a grid, we use "dch," a dict of "default choices":
- dch =  gfun.default_choices(Gr)
These are initialized in gfun.default_choices, but you typically override some of them in start_grid.py when doing the initial grid specification.  The choices are saved in a pickle file:
- pickle.dump(dch, open(Gr['gdir'] + 'choices.p', 'wb'))
You can also go back and change things in dch later using Z_edit_dch.py.

Each time you run a piece of code it makes a new grid file with the name altered to indicate what happened.  The names start as: grid_m00_r00_s00_x00.nc with the letters and numbers indicating changes to: mask, river, smoothing, or extras.

You can use plot_grid.py to look at any of the grids.  You can override the default grid by calling it in ipython like:
- run plot_grid.py -g 'sal0'

Suggested order to run the code, and which dch items are used at each step:

* start_grid.py (can take awhile for big grids; be patient; e.g. 15 min for cas6)
    analytical
    t_dir/t_list
    use_z_offset/z_offset

* make_mask.py
    z_land
    unmask_coast
    remove_islands

* edit_mask.py (to get rid of obvious issues like Lake Washington and river channels)

* carve_rivers.py

* edit_mask.py (for real this time, perhaps running many times)
	NOTE: you can use an optional command line argument -d ## to change "dval"
	the carving depth used for lines or points from its default value of 5 m.

* smooth_grid.py
    use_min_depth/min_depth

* carve_rivers.py (not needed)

* smooth_grid.py (not needed)

* make_extras.py
    min_depth (enforced for whole grid)

* grid_to_LiveOcean.py (should be run on fjord for large grids)
    nudging_edges
    nudging_days
	(also saves S-coordinate info, so be careful what file you choose in the code)

