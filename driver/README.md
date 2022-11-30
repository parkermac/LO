# README for driver

### This code is the drivers for creating forcing, running ROMS, post-processing, and other jobs that require us to control complex, repetitive tasks, either for a single forecast or for a long series of days.

---
#### `driver_forcing3.py`

This runs any of the forcing jobs, for one or more days, for any [gridname].

---

#### `driver_roms1.py` **OBSOLETE**

This runs ROMS for a single forecast or for many days. It is organized to use the LO run naming system: [gtagex] = [gridname]\_[tag]\_[ex_name].

This is much improved from the LiveOcean version:
- Being implemented in python instead of a shell script makes the code easier to follow.
- There are a number of optional flags to allow testing of each of the elements.
- The screen output is cleaned up, and there is a "verbose" mode.
- We use Path objects.
- It can run on both **klone** and **mox** in the hyak system.  You need to set specifics about each machine for each user in `LO_user/get_lo_info.py`.
- There is a new [tag_alt] logic which allows you to add new forcing variations to some existing [gtag], but do all the ROMS writing to a different [gtag]. The [tag_alt] flag refers to the existing one where the forcing is read from.
- Runs forecast as three separate days.
- Saves blowup log and last history file (and cleans these up later).
- Adds timestamp and history file number info to stdout when blowup happens.
- This is the first use of the newly-reorganized batch files, kept in `LO/driver/batch`, and with numbers (e.g klone1_) corresponding to the driver.
- This also works with the updated ROMS executables such as uu0k. It knows where to look for them by seeing if the first letter of the ex_name is repeated.
- No longer requires the -s kwarg unless the start_type is new (because the default is continuation).

---

#### `driver_roms2.py` **OBSOLETE**

Updated version of driver_roms1.py.
- fixed a bug in handling of "new" start type
- streamlined screen output
- only works with updated ROMS (LO_roms_source)

---

#### `driver_roms3.py` **USE THIS ONE**

Updated version of driver_roms2.py.
- Assumes forcing is in [gridname]
- Better error checking on mox and klone. More reliable.
- Works with perfect restart.

---

#### `driver_post1.py`

This is for running all the post-processing jobs.  It is mainly aimed at the daily forecast.  It checks to see that all expected history files are in place before beginning.

It is designed to work with a forecast that shows up as three separate days, as would be produced by `driver_roms1.py`.

---

`timestep_record.py` is a utility to make a pandas Series from the timesteps used in a model run over an arbitrary time span.

`test_Ldir.py` is a tool to test what is in Ldir. It avoids the use of the lo_tools package so it can run on the hyak clusters (mox and klone), and cron.

`test_loenv.py` also tests what is in Ldir, but using the lo_tools package.

---

#### batch

This folder contains the files used to create the batch files used to run ROMS on klone or mox.  It is meant to replace the files in LO/dot_in/shared because the batch file naming only interacts with the driver_roms code, not with the dot_in_code.  Each case consists of 2 files, named for the cluster (klone or mox) and the driver number - e.g. these files go with `driver_roms1.py`:

- `klone1_batch_BLANK.sh` is the template
- `klone1_make_batch.py` fills in the template (when called by `driver_roms1.py`)

and these produce `LO_roms/[gtagex]/[f_string]/klone1_batch.sh`

---

#### crontabs

This folder contains text files of the current crontabs for all computers in the system that use them.

You can install one with a command (from LO/driver) such as:

```
crontab crontabs/apogee.txt
```

The LANG line in the klone crontab appears to be needed when the driver starts to copy the forcing files from a remote machine like apogee.  In the past I also had to set HOSTNAME=klone, but I believe the newer drivers get around this problem. I suspect the quotes around the LOdnew entry are not needed.

Note that even though mox1 and mox2 run the same compute cluster, the crontab is specific to which one you were logged into when you created it.  So if I go to mox2 the crontab is empty.  For this reason I always use a mox1 alias when logging onto mox.
