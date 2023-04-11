# USER NOTES

## This is about blog of changes I make to the LO and related code that other user might want to be aware of.

- 2022.05.24 I added lines to `LO/dot_in/ae0_v0_uu1k/make_dot_in.py` that enable it to find lo_tools modules even when the code it run in LO_user. If you have copied this code to make your own dot_in instance in LO_user you will need to add these lines. The reason we have to add this path at all is that when running on klone we are using a generic python3 environment, not loenv, and so lo_tools is not automatically available.
- 2022.06.06 I changed the logic in `LO/driver/driver_roms2.py` for taking a nap during the forecast time period. Use git pull on klone to get these changes.
- 2022.06.08 I edited `LO/driver/driver_roms2.py` to try to catch sbatch errors that resulted in no log file being written.  Now this treats this as a blowup, and should continue. Use git pull on klone to get these changes.
- 2022.06.08 I edited `LO/driver/driver_forcing.py` to make "continuation" the default "start_type", hence unless you specifically need a "new" start type you can drop the "-s continuation" when you run this code from the command line. Use git pull on apogee or perigee to get these changes.
- 2022.06.10 I changed the times in `LO/driver/driver_roms2.py` for taking a nap during the forecast time period. Use git pull on klone to get these changes.
- 2022.06.15 I added an LO_user hook to `LO/plotting/pan_plot.py` to look for a user version of `roms_plots.py`.
- 2022.09.14 I made an edit to the klone batch script template to fix a bug introduced by the hyak maintenance yesterday. Just do `git pull` in LO to get this.

---

2022.10.18 While testing Aurora's new traps00/rivers.nc file it occurred to me that a cleaner way to handle the forcing is to collect it by [gridname] instead of [gridname]_[tag]. Then use [tag] exclusively to indicate choices made by a dot_in instance. The reason this is helpful is that you don't have to mess around with renaming the forcing collection just to be consistent with a new tag. Much cleaner!
- Result: all forcing ends up in, for example: LO_output/forcing/cas6/[frc]
- Required changes:
	- forcing code:
		- edit all the "00" code to reflect the new output location.
		- Also create lo_tools/forcing_argfun2.py so that the changes don't interfere with the current forecast. Note that the forcing code will no longer accept a tag argument.
		- Also create driver/driver_forcing3.py to use the new system
		- To do: implement similar changes for the Analytical and Nesting forcing and dot_in code.
	- driver_roms3.py: change where it looks for (and where it puts) the forcing files, and remove tag_alt machinery. I also cleaned up the sbatch code, making "3" versions for both klone and mox.
	- dot_in: change path to forcing, starting with cas6_v00_uu0mb
	- On apogee, move LO_output/forcing/cas6_v00 to LO_output/forcing/cas6

---

2023.03.14 I have done a lot of work on processing of observational cast data (bottle and ctd). You can see current notes in the README for LO/obs, and be sure to update your processed files using those currently in LO_output/obs on perigee or apogee. I have not yet updates the files in LO_data/obs that were orignally used to create these processed products. Currently the "official versions" are on my laptop.

---

2023.04.11 Over the weekend I started a new run on mox, and along the way made a change to how I organize the bio code that may be useful. In LO_roms_user/x1b I added the fennel.h code to the x1b folder, instead of leaving it in LO_roms_source_alt. Then I changed a line in build_roms.sh so that MY_ANALYTICAL_DIR points to x1b as well. The motivation is that it keeps any edits of fennel.h connected to a specific [ex_name].
