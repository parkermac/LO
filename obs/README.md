# README for obs

#### These folders contain the highly specific code for processing observations into LO-standardized formats that can be used for model validation and other analysis tasks.

---

## Organization

#### Naming Notes
- [source] is some string identifier of the database, e.g. 'dfo', which is then further identified by
- [otype] observation type, e.g. 'ctd' or 'bottle', and finally by
- [year] e.g. 2017.


#### Input/Output file naming conventions
- Input: LO_data/obs/[source]/.csv, .xls, etc.
- Code: LO/obs/[source]/process_ctd_bottle.py or whatever is called for
- Output: LO_output/obs/[source]/[otype]/[year].p and info_[year].p or whatever is called for

---

## CTD and Bottle output conventions

#### info_[year].p is a pickled DataFrame that has (x,y,t) for each cast in that year, with

**index** = cid: a unique number or string that identifies a single cast in the corresponding ctd or bottle DataFrame

**columns**:
- lon [degrees -180:180]
- lat [degrees -90:90]
- time [datetime UTC]
- name [station name or None]

#### [year].p is a pickeld DataFrame with all the processed data for all the casts.

**index** is just numbers, one for each depth

**columns**:
- 'cid', 'lon', 'lat', 'time', name are IDENTICAL for all rows of a given cast, and correspond to a single row in the info_[year].p DataFrame
- 'z' [z-position in meters, positive up, zero at free surface]
- 'salt (SA g kg-1)' [Absolute Salinity in g kg-1]
- 'temp (CT degC)' [Conservative Temperature in degC]
- 'DO (uM)' [dissolved oxygen in umol L-1]
- and whatever other data columns there are...

Note that we use a sort of cumbersome naming convention for the data columns, including the units explicitly so it is clear what we are dealing with in subsequent use. After I have done more processing I will create a master list of preferred and acceptable variable names, e.g. we might allow for both 'salt (SA g kg-1)' and 'salt (psu)'.

---

## Mooring output conventions

_Under construction_

---

#### Sources

- **dfo** Canadian ctd and bottle data
- _other sources have not yet been brought up to current organizational standards..._
