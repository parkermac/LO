# Notes on the files for external users

#### `driver_post1.py` runs a sequence of post-processing jobs, and several of these are made specifically for external users.

In all cases the output files are copied to the APL server [https://liveocean.apl.uw.edu/output/](https://liveocean.apl.uw.edu/output/) every morning starting around 6 AM Pacific time (if all goes according to plan) where they are publicly available.  Each day's files are in a folder with a name like f2022.01.11.

The lineup of extraction jobs can be found in `LO/driver/driver_post1.py`.  For example:
```
job_list = ['nest_wgh', 'surface1', 'layers1', 'ubc1', 'sequim1', 'critfc1',
    'daymovie0', 'drifters0', 'archive0']
```
and you can look in LO/post/[job name] to see how the processing is handled.  The number in the job name refers to the version, and it may change as I update things, but on the APL server the output file name omits that number, and should remain stable.  For example, surface[whatever number] will always make `surface.nc`.

There is also a small text file [jobname]_done.txt generated and sent to the server when each file is finished being written.  For example, on the day I wrote this, the surface_done.txt file contains the string "2022.01.11 05:47:34" (local time when it was written).  Users could look for this file before downloading their extraction.

#### Many of the jobs have people who are expecting the results. Below is the lineup, the users if any, and other info.

---

### nest_wgh

This generates the ocean forcing expected by our nested sub-model of Willapa Bay and Grays Harbor.

---

### surface1

This is an extraction of surface fields served by IOOS.

Output: `surface.nc` [257 MB]

Contacts:
- Tim Giguere <Tim.Giguere@rpsgroup.com> IOOS EDS Viewer
- Patrick Tripp <Patrick.Tripp@rpsgroup.com> IOOS S3-THREDDS server
- Kelly Knee <Kelly.Knee@rpsgroup.com> Admin.

---

### layers1 and layers2

This is an extraction of fields (including processed ones like Aragonite Saturation State) on the surface, the bottom, and several horizontal depth layers. It is what appears on the NANOOS NVS, and is used by an aquaculture services company, SCOOT.

Output: `layers.nc` [1.7 GB]

Contacts at NANOOS:
- Alex Dioso <adioso@uw.edu>
- Troy Tanner <troyt@apl.washington.edu>
- Roxanne J Carini <rjcarini@uw.edu>
- Craig Risien <craig.risien@oregonstate.edu>
- Jan A. Newton <janewton@uw.edu> Admin.

Contacts at SCOOT:
- Evan Goodwin <evan@scootscience.com>
- Iwen Su <iwen@scootscience.com>
- Connor Dibble <connor.dibble@scootscience.com>

---

### layers_uv

This is an extraction of velocity on layers in JdF for some glider operators at APL. The results go to an AWS s3 bucked on kopah.

Contact:
- Peter Brodsky <pmb@uw.edu>

---

### ubc

This is a small extraction of model fields at the mouth of Juan de Fuca, used to make the open boundary condition for Susan Allen's Salish SeaCast model.

Output: `ubc.nc` [10 MB]

Contacts:
- Susan Allen <sallen@eoas.ubc.ca>
- Doug Latornell <djlatornell@gmail.com>

---

### sequim

This is a small extraction of model fields near Sequim Bay, used to make the open boundary condition for a model run by Zhaoqing Zhang at PNNL.

Output: `sequim.nc` [54 MB]

Contacts:
- Zhaoqing Yang <Zhaoqing.Yang@pnnl.gov>,
- Clair Ma <clair.ma@pnnl.gov>

---

### harcourt [not currently running 2023.05.18]

This is an extraction of model fields off the WA coast for an experiment on sound speed being conducted by Ramsey Harcourt and colleagues at APL.

Output: `harcourt.nc` [1 GB] and a folder `harcourt_hourly`

Contacts:
- Ramsey Harcourt: <harcourt@uw.edu>

---

### critfc [not running as of 2024.08.19]

This is a large extraction of model fields used by Charles Seaton as backup boundary conditions for his CRITFC model of the Columbia River.

Output: `critfc_[salt,temp].nc` (two files) [2.8 GB total]

Contact:
- Charles Seaton <cseaton@critfc.org>

---

### daymovie

This creates the daily forecast movies that end up in Parker's LiveOcean website. The most critical one is `Phab_full_salt_top.mp4` which is used by the people that make decisions about the coastal razor clam harvest.

Contacts:
- Ryan McCabe <ryan.mccabe@noaa.gov>
- Dan Ayres <Daniel.Ayres@dfw.wa.gov>
- Matt Hunter <Matthew.V.HUNTER@odfw.oregon.gov>

---

### drifters

This runs the particle tracking runs using the tracker tool, and then sends them a JSON to Parker's website for the interactive javascript tools there.

---

### archive [not currently running 2024.08.25]

This copies the forecast history files to a standard place on perigee. Some users that have an account there use them, e.g. Harper Simmons who uses them to create boundary conditions for a nested forecast of Hood Canal.

Output: on perigee in /data1/parker/LO_roms/cas6_v1_live

Contact:
- Harper Simmons <harpers@uw.edu>

---

### azure [not currently running 2023.05.18]

This copies the forecast history files to blob storage on azure. Is is a custom job for Rich Signell to access during a coding event.

Contact:
- Rich Signell <rsignell@usgs.gov>

---
