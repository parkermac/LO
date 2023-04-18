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
- 'DO (uM)' [dissolved oxygen in umol L-1 = mmol m-3 = "micro-molar"]
- We also use units of (uM) for NO3, NO2, NH4, PO4, SiO4, TA, and DIC
- 'Chl (mg m-3)' [Chlorophyll a]

Note that we use a sort of cumbersome naming convention for the data columns, including the units explicitly so it is clear what we are dealing with in subsequent use. After we have done more processing I will create a master list of preferred and acceptable variable names.

#### Sources that have been processed up to these specifications

Note: After you download these to your laptop you can get a quick look at the data by running LO/obs/plot_ctd_bottle.py.

Folders on perigee in `/data1/parker/LO_output/obs` and on apogee in `/dat1/parker/LO_output/obs`:

- **dfo1** Canadian data from a large database created by Susan Allen's group. The original data came from: https://data.cioospacific.ca/erddap/index.html. This replaces the older "dfo" that was parsed from an SQL database.
  - bottle 1930-2021 with gaps in the early years, has some Chl data
  - ctd 1965-2021
- **ecology** Department of Ecology monthly repeat station data. Covers 2008-2017 for bottles and 2008-2019 for ctd casts. 300-400 casts per year, monthly at 39 stations in Puget Sound and the coastal estuaries. Processing note: I added SA, CT, and DO (uM) to the processed bottle files using values interpolated from ctd casts.
  - bottle 2008-2017, added SA, CT, DO from ctd casts
  - ctd 2008-2019, has some Chl data
- **nceiSalish** Data from WOAC and other cruises, mainly in Puget sound and JdF. The original data came from: https://www.ncei.noaa.gov/access/ocean-carbon-acidification-data-system/oceans/SalishCruise_DataPackage.html.
  - bottle 2008-2018, good DIC and TA, 40-199 casts per year
  - ctd None
- **nceiCoastal** Data from WCOA and other North American coastal cruises with carbon data. Good for carbon data on the WA/OR shelf. The original data came from: https://www.ncei.noaa.gov/data/oceans/ncei/ocads/metadata/0219960.html.
  - bottle 2011-2017 in the LO model domain, good DIC and TA, 41-56 casts per year
  - ctd None

Notes on usage:
- 2017 is a good year for validation, with bottle coverage from all four sources.
- If you are interested in Chl, you will only find it in dfo1-bottle and ecology-ctd.
- For ctd data (much finer vertical resolution, and potentially more info about bottom hypoxia) there are only two sources: dfo1 and ecology.
- For carbon data (DIC and TA) you will only find them in nceiSalish and nceiCoastal, hence no coverage in Canadian waters.

---

## Mooring output conventions

_Under construction_
