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

---

#### `driver_post.py`

This is for running all the post-processing jobs.  It is mainly aimed at the daily forecast.  It checks to see that all expected history files are in place before beginning.

---

`copy_to_azure.py` is a utility to copy a file like a box extraction to azure for sharing.

`test_Ldir.py` is a tool to test what is in Ldir. It avoids the use of the lo_tools package so it can run on the hyak clusters (mox and klone), and cron.
