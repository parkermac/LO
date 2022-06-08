# USER NOTES

## This is about blog of changes I make to the LO and related code that other user might want to be aware of.

- 2022.05.24 I added lines to `LO/dot_in/ae0_v0_uu1k/make_dot_in.py` that enable it to find lo_tools modules even when the code it run in LO_user. If you have copied this code to make your own dot_in instance in LO_user you will need to add these lines. The reason we have to add this path at all is that when running on klone we are using a generic python3 environment, not loenv, and so lo_tools is not automatically available.
- 2022.06.06 I changed the logic in `LO/driver/driver_roms2.py` for taking a nap during the forecast time period. Use git pull on klone to get these changes.
- 2022.06.07 I edited `LO/driver/driver_roms2.py` to try to catch sbatch errors that resulted in no log file being written.  Now this treats this as a blowup, and should continue. Use git pull on klone to get these changes.
