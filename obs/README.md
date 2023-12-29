# README for obs

#### These folders contain the highly specific code for processing observations into LO-standardized formats that can be used for model validation and other analysis tasks.

---

## Organization

#### Naming Notes
- [source] is some string identifier of the database, e.g. 'dfo', which is then further identified by
- [otype] observation type, e.g. 'ctd', 'bottle', or 'moor', and finally by
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

Here is an example of the contents of one of the Ecology DataFrames:

```
cid       time        lat         lon    name  ...  NO2 (uM)  NH4 (uM)  PO4 (uM)  SiO4 (uM)  cruise
0      0.0 2017-01-11  47.092040 -122.918197  BUD005  ...  0.388398  0.668405  2.632084  84.497398    None
1      0.0 2017-01-11  47.092040 -122.918197  BUD005  ...  0.340418  0.726122  2.629354  82.600006    None
2      0.0 2017-01-11  47.092040 -122.918197  BUD005  ...  0.348933  0.628245  2.605122  82.263443    None
3      1.0 2017-01-11  47.276485 -122.709575  CRR001  ...  0.103605  0.037389  2.551934  71.561836    None
4      1.0 2017-01-11  47.276485 -122.709575  CRR001  ...  0.117088  0.051287  2.561856  70.790131    None
..     ...        ...        ...         ...     ...  ...       ...       ...       ...        ...     ...
939  276.0 2017-12-13  47.398148 -122.929592  HCB007  ...  0.303905  1.652716  3.884351  84.559181    None
940  277.0 2017-12-13  47.670000 -122.820000  HCB010  ...  0.318795  0.071606  1.763635  77.125214    None
941  277.0 2017-12-13  47.670000 -122.820000  HCB010  ...  0.090017  0.102635  2.565202  61.486889    None
942  277.0 2017-12-13  47.670000 -122.820000  HCB010  ...  0.042560  0.000384  2.714923  58.599747    None
943  277.0 2017-12-13  47.670000 -122.820000  HCB010  ...  0.120723  0.001253  2.887872  57.303364    None
```

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
- **collias** Historic Data collected by Eugene Collias. 368,858 results from the years 1932 until 1975. Data was collected using approximately 15-20 different methods and captures 9 different parameters. Data came from https://apps.ecology.wa.gov/eim/search/Eim/EIMSearchResults.aspx?ResultType=EIMStudyTab&LocationWRIAs=2 (Select the dropdown menu for Collias and then select Download. This tirggers an EIM request and eventually they get it to you.) And here is a source of more information: https://apps.ecology.wa.gov/eim/search/Detail/Detail.aspx?DetailType=Study&SystemProjectId=99971885 with links to some great studies.
  - bottle 1932-1975. Missing 1943-1948 and 1973. Good SA, CT and DO, then less for NO3, NO2, and SiO4.
  - ctd None.

Notes on usage:
- 2017 is a good year for validation, with bottle coverage from all four sources.
- If you are interested in Chl, you will only find it in dfo1-bottle and ecology-ctd.
- For ctd data (much finer vertical resolution, and potentially more info about bottom hypoxia) there are only two sources: dfo1 and ecology.
- For carbon data (DIC and TA) you will only find them in nceiSalish and nceiCoastal, hence no coverage in Canadian waters.

---

## Mooring output conventions

#### Input/Output file naming conventions
- Input: LO_data/obs/[source]/[whatever]
- Code: LO/obs/[source]/process_data.py or whatever is called for
- Output: LO_output/obs/[source]/[otype]/[whatever name makes sense].nc

#### Output format
We store processed mooring files in NetCDF, created from xarray Datasets.
- Each variable is an array packed with dimensions (time,z).
- the time coordinate is best made as a pandas DatetimeIndex object, with regular interval.
- z is a vector of vertical position, packed bottom to top
- The variable naming convention follows that of the bottle and ctd data. We also store metadata about each variable, which can be accessed for example as ds.SA.attrs['long_name'] or ds.SA.attrs['units']

Here is an example of the contents of processed ORCA data `LO_output/obs/orca/moor/DB_daily.nc`:
```
<xarray.Dataset>
Dimensions:   (time: 4243, z: 104)
Coordinates:
  * time      (time) datetime64[ns] 2010-06-08 2010-06-09 ... 2022-01-18
  * z         (z) float64 -105.1 -104.1 -103.1 -102.1 ... -4.958 -3.966 -2.975
Data variables:
    SA        (time, z) float64 nan nan nan nan nan ... 26.79 25.35 23.15 20.4
    CT        (time, z) float64 nan nan nan nan nan ... 7.696 7.43 7.411 7.204
    DO (uM)   (time, z) float64 nan nan nan nan nan ... 267.6 289.5 308.5 317.5
    NO3 (uM)  (time, z) float64 ...
    FLUOR     (time, z) float64 ...
    PAR       (time, z) float64 ...
    SIG0      (time, z) float64 ...
Attributes:
    Station Name:  Dabob Bay
    lon:           -122.8029
    lat:           47.8034
```
#### Sources that have been processed up to these specifications

- **orca** ORCA profiling mooring data from NANOOS. Six moorings around Puget Sound, processed into regular daily profiles (with gaps) by Erin Broatch. See Erin's report in LO_data/obs/ORCA/orca_report.pdf for the original source at NANOOS, however the link to that source is broken, but this should give the same access: https://nwem.apl.washington.edu/prod_DataReq.shtml#. Now (as of late 2023) there is also a new automated server at https://data.nanoos.org/erddap/tabledap/index.html?page=1&itemsPerPage=1000. Eventually we will want to do the processing from that data source.
 - time range 2005 to the end of 2021 (at most), daily, with many gaps
 - variables: SA, CT, DO (uM), and SIG0 are pretty robust. NO3 (uM), FLUOR, and PAR are carried along in the processing, but I am not sure they can be trusted, or that the even have much data.
 - There is more plotting code to help explore the data in different ways in LO/obs/orca.
![ORCA Mooring Data](readme_plots/orca.png)
