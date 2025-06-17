# README for pre/glorys

### The code in this folder is designed to allow you to download GLORYS ocean fields which are then used to create forcing files for ROMS backfill and forecast runs. Development notes are in LPM/glorys.

General info about glorys products:

https://data.marine.copernicus.eu/products

Forecast info:

https://data.marine.copernicus.eu/product/GLOBAL_ANALYSISFORECAST_PHY_001_024/description

Click on: Data Access -> Dataset -> Daily -> Form to find specific variables, units, and ranges.

Hindcast info:

https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_PHY_001_030/description

Data info:
- 50 vertical layers
- regular 8 km (1/12 degree) grid
- daily = daily average centered at midnight at start of that day

Daily fields are available over three different time ranges, accessed using different product_id strings. I call these **products** 'hindcast', 'interim', and 'forecast' in my code. For hindcast and interim you get all physical variables in one file. For forecast you get each separately.

---

**lo_tools/glorys_functions.py** is a module of useful functions for downloading and processing (e.g. interpolate to ROMS grid) glorys fields. An important function is **get_glorys_region()** which returns a list of lon-lat limits associated with a **region** name such as "region7" for a region that encompasses the cas7 domain.

---

**get_glorys_days.py** is the main program to get the glorys data. 

- Input: command line arguments for start day, end day, product, and region.

- Output for hindcast and interim products:
```
LO_data/glorys/[region]/hindcast_[date string].nc
```

- Output for forecast product:
```
LO_data/glorys/[region]/forecast_[date_string]/
  forecast_[so,thetao,cur,zso].nc
```
Note that "cur" contains velocity variables uo and vo.

---

**get_time_ranges.py** will print out the current time ranges available for each product. As of 2025.06.16 these were:

```
hindcast
cmems_mod_glo_phy_my_0.083deg_P1D-m
INFO - 2025-06-16T15:22:57Z - Selected dataset version: "202311"
INFO - 2025-06-16T15:22:57Z - Selected dataset part: "default"
Range: 1993.01.01 to 2021.06.30

interim
cmems_mod_glo_phy_myint_0.083deg_P1D-m
INFO - 2025-06-16T15:23:14Z - Selected dataset version: "202311"
INFO - 2025-06-16T15:23:14Z - Selected dataset part: "default"
Range: 2021.07.01 to 2025.05.27

forecast
cmems_mod_glo_phy-so_anfc_0.083deg_P1D-m
INFO - 2025-06-16T15:23:28Z - Selected dataset version: "202406"
INFO - 2025-06-16T15:23:28Z - Selected dataset part: "default"
Range: 2022.06.01 to 2025.06.22
```