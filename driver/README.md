# README for driver

### This code is the drivers for creating forcing, running ROMS, post-processing, and other jobs that require us to control complex, repetitive tasks, either for a single forecast or for a long series of days.

---

#### `driver_roms3.py`

This is the primary driver for running ROMS on both klone and mox. It can run both daily forecasts and long backfill runs. Read through the comments at the top of the code for more details and development notes. It has many loops built in to try an operation multiple times to make it more reliable.

---

#### `driver_romsN.py`

This is nearly identical to driver_roms3.py but is customized to do nested runs. The differences are so minor I think it could be absorbed into driver_roms3.py without difficulty.

---

#### `driver_post1.py`

This is for running all the post-processing jobs.  It is mainly aimed at the daily forecast.  It checks to see that all expected history files are in place before beginning.

It is designed to work with a forecast that shows up as three separate days, the current standard.
---

#### `driver_post2.py`

This is for running more post-processing jobs.  It is mainly aimed at the daily forecast.  It checks to see that all expected history files are in place before beginning.

This was created to have a different lineup of jobs, aimed at the nested forecast for Willapa Bay and Grays Harbor. Otherwise identical to `driver_post1.py`.

---

#### `driver_post_backfill.py`

This runs post-processing jobs (currently only layers1) as multi-day backfill jobs. It only exists as a convenience for generating extractions that some users wanted.

---

#### Miscellany

`timestep_record.py` is a utility to make a pandas Series from the timesteps used in a model run over an arbitrary time span.

`test_Ldir.py` is a tool to test what is in Ldir. It avoids the use of the lo_tools package so it can run on the hyak clusters (mox and klone), and cron.

`test_loenv.py` also tests what is in Ldir, but using the lo_tools package.

---

#### batch

This folder contains the files used to create the batch files used to run ROMS on klone or mox.  It is meant to replace the files in LO/dot_in/shared because the batch file naming only interacts with the driver_roms code, not with the dot_in_code.  Each case consists of 2 files, named for the cluster (klone or mox) and the driver number - e.g. these files go with `driver_roms3.py`:

- `klone3_batch_BLANK.sh` is the template
- `klone3_make_batch.py` fills in the template (when called by `driver_roms3.py`)

and these produce `LO_roms/[gtagex]/[f_string]/klone_batch.sh`

---

#### crontabs

This folder contains text files of the current crontabs for all computers in the system that use them.

You can install one with a command (from LO/driver) such as:

```
crontab crontabs/apogee.txt
```

The LANG line in the klone crontab appears to be needed when the driver starts to copy the forcing files from a remote machine like apogee.  In the past I also had to set HOSTNAME=klone, but I believe the newer drivers get around this problem. I suspect the quotes around the LOdnew entry are not needed.

Note that even though mox1 and mox2 run the same compute cluster, the crontab is specific to which one you were logged into when you created it.  So if I go to mox2 the crontab is empty.  For this reason I always use a mox1 alias when logging onto mox.
