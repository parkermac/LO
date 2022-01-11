# Notes on the files for external users

#### `driver_post1.py` runs a sequence of post-processing jobs, and several of these are made specifically for external users.

In all cases the output files are copied to the APL server [https://liveocean.apl.uw.edu/output/](https://liveocean.apl.uw.edu/output/) every morning starting around 6 AM Pacific time (if all goes according to plan) where they are publicly available.  Each day's files are in a folder with a name like f2022.01.11.

The lineup of extraction jobs can be found in LO/driver/driver_post1.py.  Currently it is:

- 'surface1', 'layers1', 'ubc1', 'sequim1', 'critfc1'

and you can look in LO/post/[job name] to see how the processing is handled.  The number in the job name refers to the version, and it may change as I update things, but on the APL server the output file name omits that number, and should remain stable.  For example, surface[whatever number] will always make `surface.nc`.

There is also a small text file [jobname]_done. generated when each file is finished being written.  For example, on the day I wrote this, the surface_done.txt file contains the string "2022.01.11 05:47:34" (local time when it was written).  Users could look for this file before downloading their extraction.

---

### surface

Output: `surface.nc` [257 MB]

Contacts:
- Tim Giguere <Tim.Giguere@rpsgroup.com> IOOS EDS Viewer
- Patrick Tripp <Patrick.Tripp@rpsgroup.com> IOOS S3-THREDDS server
- Kelly Knee <Kelly.Knee@rpsgroup.com> Admin.

---

### layers

Output: `layers.nc` [1.7 GB]

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

### ubc

Output: `ubc.nc` [10 MB]

Contacts:
- Susan Allen <sallen@eoas.ubc.ca>
- Doug Latornell <djlatornell@gmail.com>

---

### sequim

Output: `sequim.nc` [54 MB]

Contacts:
- Zhaoqing Yang <Zhaoqing.Yang@pnnl.gov>,
- Clair Ma <clair.ma@pnnl.gov>

---

### critfc

Output: `critfc_[salt,temp].nc` (two files) [2.8 GB total]

Contact:
- Charles Seaton <cseaton@critfc.org>
