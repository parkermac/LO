# Notes on the files for external users

#### `driver_post.py` runs a sequence of post-processing jobs, and several of these are made specifically for external users.

In all cases the output files are copied to boiler, where all groups have accounts so they can scp.  The output folder is:

(*) = /data1/parker/LiveOcean_roms/output/cas6_v3_lo8b/f[date string]

After each file is copied there a small text file is generated so that users can know it is ready.  It is of the form:

(*)/[surface,layers,ubc,critfc]_done.txt

and contains a text string of the form "2021.11.11 08:29:00" (local time when it was written).

---

### surface0

Output: (*)/ocean_surface.nc

Contacts:
- Tim Giguere <Tim.Giguere@rpsgroup.com>
- Patrick Tripp <Patrick.Tripp@rpsgroup.com>
- Kelly Knee <Kelly.Knee@rpsgroup.com>

---

### layers0

Output: (*)/ocean_surface.nc

Contacts:
- Susan Allen
- Doug Latournell

---

### ubc0

Output: (*)/low_passed_UBC.nc

Contacts:
- Susan Allen
- Doug Latournell

---

### critfc0

Output: (*)/cmop_[salt,temp]_nu.[YYYY-MM-DD].nc (two files)

Contacts:
- Charles Seaton
