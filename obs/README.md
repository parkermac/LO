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
- cruise [cruise name if available]

#### [year].p is a pickled DataFrame with all the processed data for all the casts.

**index** is just numbers, one for each depth

**columns**:
- 'cid', 'lon', 'lat', 'time', 'name', 'cruise' are IDENTICAL for all rows of a given cast, and correspond to a single row in the info_[year].p DataFrame
- 'z' [z-position in meters, positive up, zero at free surface]
- 'SA' [Absolute Salinity in g kg-1]
- 'CT' [Conservative Temperature in degC]
- 'DO (uM)' [dissolved oxygen in umol L-1]
- and whatever other data columns there are,

Note that we use a sort of cumbersome naming convention for the data columns, including the units explicitly so it is clear what we are dealing with in subsequent use. After we have done more processing I will create a master list of preferred and acceptable variable names.

---

## Mooring output conventions

_Under construction_

---

#### Sources that have been processed up to these specifications

Note: After you download these to your laptop you can get a quick look at the data by running LO/obs/plot_ctd_bottle.py.

Folders on perigee in `/data1/parker/LO_output/obs`)

- **dfo** Canadian data from a large SQL database created by Susan Allen's group. Covers 1930-2019, with some gaps. Currently only bottle data has been processed.
- **ecology** Department of Ecology monthly repeat station data. Covers 2008-2017 for bottles and 2008-2019 for ctd casts. 300-400 casts per year, monthly at 39 stations in Puget Sound and the coastal estuaries. Processing note: I added SA, CT, and DO (uM) to the processed bottle files using values interpolated from ctd casts.
- **nceiSalish** Data from WOAC and other cruises, mainly in Puget sound and JdF, 2008-2018. 40-100 casts per year. The original data came from: https://www.ncei.noaa.gov/access/ocean-carbon-acidification-data-system/oceans/SalishCruise_DataPackage.html. Bottles only.
- **nceiCoastal** Data from WCOA and other North American coastal cruises with carbon data. Covers 2011-2017 in the LO model domain, with some gaps, 41-56 casts per year. Good for carbon data on the WA/OR shelf. The original data came from: https://www.ncei.noaa.gov/data/oceans/ncei/ocads/metadata/0219960.html. Bottles only.
