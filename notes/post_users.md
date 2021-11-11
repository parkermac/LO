# Notes on the files for external users

#### `driver_post.py` runs a sequence of post-processing jobs, and several of these are made specifically for external users.

In all cases the output files are copied to boiler, where all groups have accounts so they can get them via scp.  For historical consistency, the output folder is:

(*) = `/data1/parker/LiveOcean_roms/output/cas6_v3_lo8b/f[date string]`

even though the actual run [gtagex] has changed.  After each file is copied there a small text file is generated so that users can know it is ready.  It is of the form:

`(*)/[surface,layers,ubc,critfc]_done.txt`

containing a text string of the form "2021.11.11 08:29:00" (local time when it was written).

---

### surface0

Output: `(*)/ocean_surface.nc` (also sent to Azure) [0.2 GB]

Contacts:
- Tim Giguere <Tim.Giguere@rpsgroup.com> IOOS EDS Viewer
- Patrick Tripp <Patrick.Tripp@rpsgroup.com> IOOS S3-THREDDS server
- Kelly Knee <Kelly.Knee@rpsgroup.com> Admin.

---

### layers0

Output: `(*)/ocean_surface.nc` [1.8 GB]

Contacts at NANOOS:
- Alex Dioso <adioso@uw.edu>
- Troy Tanner <troyt@apl.washington.edu>
- Roxanne J Carini <rjcarini@uw.edu>
- Craig Risien <craig.risien@oregonstate.edu>
- Jan A. Newton <janewton@uw.edu> Admin.

Contacts at SCOOT:
- Henry Amestoy <henry.amestoy@scootscience.com>
- Evan Goodwin <evan@scootscience.com>
- Iwen Su <iwen@scootscience.com>

---

### ubc0

Output: `(*)/low_passed_UBC.nc` [small]

Contacts:
- Susan Allen <sallen@eoas.ubc.ca>
- Doug Latornell <djlatornell@gmail.com>

---

### critfc0

Output: `(*)/cmop_[salt,temp]_nu.[YYYY-MM-DD].nc` (two files) [3 GB total]

Contact:
- Charles Seaton <cseaton@critfc.org>
