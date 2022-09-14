# USER NOTES

## This is about blog of changes I make to the LO and related code that other user might want to be aware of.

- 2022.05.24 I added lines to `LO/dot_in/ae0_v0_uu1k/make_dot_in.py` that enable it to find lo_tools modules even when the code it run in LO_user. If you have copied this code to make your own dot_in instance in LO_user you will need to add these lines. The reason we have to add this path at all is that when running on klone we are using a generic python3 environment, not loenv, and so lo_tools is not automatically available.
- 2022.06.06 I changed the logic in `LO/driver/driver_roms2.py` for taking a nap during the forecast time period. Use git pull on klone to get these changes.
- 2022.06.08 I edited `LO/driver/driver_roms2.py` to try to catch sbatch errors that resulted in no log file being written.  Now this treats this as a blowup, and should continue. Use git pull on klone to get these changes.
- 2022.06.08 I edited `LO/driver/driver_forcing.py` to make "continuation" the default "start_type", hence unless you specifically need a "new" start type you can drop the "-s continuation" when you run this code from the command line. Use git pull on apogee or perigee to get these changes.
- 2022.06.10 I changed the times in `LO/driver/driver_roms2.py` for taking a nap during the forecast time period. Use git pull on klone to get these changes.
- 2022.06.15 I added an LO_user hook to `LO/plotting/pan_plot.py` to look for a user version of `roms_plots.py`.
- 2022.09.14 I made an edit to the klone batch script template to fix a bug introduced by the hyak maintenance yesterday. Just do `git pull` in LO to get this.
