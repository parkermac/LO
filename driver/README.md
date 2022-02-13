# README for driver

### This code is the drivers for creating forcing, running ROMS, post-processing, and other jobs that require us to control complex, repetitive tasks, either for a single forecast or for a long series of days.

---
#### `driver_forcing.py`

This runs any of the forcing jobs, for one or more days, for any [gtag].

---

#### `driver_roms.py`

This runs ROMS for a single forecast or for many days. It is organized to use the LO run naming system: [gtagex] = [gridname]\_[tag]\_[ex_name].

This is much improved from the LiveOcean version:
- Being implemented in python instead of a shell script makes the code easier to follow.
- There are a number of optional flags to allow testing of each of the elements.
- The screen output is cleaned up, and there is a "verbose" mode.
- We use Path objects.
- It should run on both **klone** and **mox** in the hyak system.  You need to set specifics about each machine for each user in `LO_user/get_lo_info.py`.
- There is a new [tag_alt] logic which allows you to add new forcing variations to some existing [gtag], but do all the ROMS writing to a different [gtag]. The [tag_alt] flag refers to the existing one where the forcing is read from.

#### `driver_roms1.py`

This is meant to be an update to the current forecast driver `driver_roms.py`, with several improvements:

- Runs forecast as three separate days.
- Saves blowup log and last history file (and cleans these up later).
- Adds timestamp and history file number info to stdout when blowup happens.
- This is the first use of the newly-reorganized batch files, kept in `LO/driver/batch`, and with numbers (e.g klone1_) corresponding to the driver.
- This also works with the updated ROMS executables such as uu0k. It knows where to look for them by seeing if the first letter of the ex_name is repeated.
- No longer requires the -s kwarg unless the start_type is new (because the default is continuation).

---

#### `driver_post.py`

This is for running all the post-processing jobs.  It is mainly aimed at the daily forecast.  It checks to see that all expected history files are in place before beginning.

---

#### `driver_post1.py`

This is for running all the post-processing jobs.  It is mainly aimed at the daily forecast.  It checks to see that all expected history files are in place before beginning.

It is designed to work with a forecast that shows up as three separate days, as would be produced by `driver_roms1.py`.

---

`timestep_record.py` is a utility to make a pandas Series from the timesteps used in a model run over an arbitrary time span.

`test_Ldir.py` is a tool to test what is in Ldir. It avoids the use of the lo_tools package so it can run on the hyak clusters (mox and klone), and cron.

---

#### batch

This folder contains the files used to create the batch files used to run ROMS on klone or mox.  It is meant to replace the files in LO/dot_in/shared because the batch file naming only interacts with the driver_roms code, not with the dot_in_code.  Each case consists of 2 files, named for the cluster (klone or mox) and the driver number - e.g. these files go with `driver_roms1.py`:

- `klone1_batch_BLANK.sh` is the template
- `klone1_make_batch.py` fills in the template (when called by `driver_roms1.py`)

and these produce `LO_roms/[gtagex]/[f_string]/klone1_batch.sh`

---

#### crontabs

This folder contains text files of the current crontabs for all computers in the system that use them.

I think that you should be able to install one with a command (from LO/driver) such as:

```
crontab crontabs/apogee.txt
```

but I haven't tested it yet.
