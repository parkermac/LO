# README for the pre/river1 code

### This code is for gathering and processing river data, and one of the uses for the data is forcing ROMS hindcasts. This replaces pre/river.

#### It relies on functions in the module `LO/lo_tools/lo_tools/river_functions.py`, which go out and get data from external sources: USGS, NWS Forecasts, and Environment Canada websites.  The reason for putting these in lo_tools is that they are also used by the forecast system.

#### The folder was started 2023.04.09 based on pre/river. The main improvement is to simplify the process of updating the historical and climatological data every year. We also making the naming more generic, so that there are not hard-coded date ranges in the river forcing code, for example. Any river forcing code will have to be modified to make use of this new naming.

---

#### CURRENT OUTPUT STATE 2023.04.09
Collection lo_base has:
- historical flow and  flow climatology for 1980-2022
- historical temperature and temperature climatology for 1980-2020
- On 2023.05.01 I used the new pgrid/create_river_tracks.py tool to define new tracks for the four rivers entering Willapa Bay and Grays Harbor: naselle, willapa, chehalis, and humptulips. These are saved in the collection "lo_new". I then saved the originals in lo_base with "_ORIG" appended and copied the new tracks into lo_base. These are just simpler tracks that will allow a better wgh0 grid.

---

The output goes into:

(*) = LO_output/pre/river1/[ctag]

where [ctag] is a "collection tag".

---

`make_river_info.py` organizes river info from different sources.

INPUT:

Information from Sarah Gidding, Ecology, and other places.

OUTPUT:

- (*)/river_info.csv - names, gauges, scaling factors
- (*)/river_info.p - pickled DataFrame version of river_info.csv
- (*)/tracks/[river name].p - lon, lat tracks for carving river channels in a grid

NOTE: Currently this is hard-coded to the reflect the collection of rivers in LO as of 2023.04.09 (without the TRAPS rivers - they are added later). This has **ctag = lo_base**.

---

`make_historical.py` - the main tool for getting historical flow data. It accepts command line arguments for the ctag and year range. The default behavior **with no arguments** is to add the most recent year (e.g. 2022 if run in 2023) to the existing ALL_flow.p, moving the old one to ALL_flow_prev.p. This is really useful because it means you can update the historical archive with one command, especially if you do it early in the year (say January-April) before we drift into the dreaded EC Data Gap. If you do it late in the year you may have to update the ROMS river forcing extraction and edit the code to use it.

INPUT:

Raw data from XML (USGS), or scraped from html (EC), using the functions in river_functions.py.  For the USGS data you can simply hand it a date range and it will quickly return 40 years of daily flow data.  For EC data it is more complicated.  Typically over a range in the last 18 months you can get data using a date range, just like for USGS.  Earlier than this you need to scrape the "historical" data from a table on a web page, a year at a time, but these only exist up to two years ago, e.g. 2020 when it is 2022. If you do this after mid-year there is a gap in the EC data, e.g. doing this in November 2022 there is a gap for the early months of 2021. A workaround is to use as-run river extractions from LiveOcean forcing files, which are created using code in LO/extract/river.

OUTPUT:

(++) = (*)/Data_historical/ALL_flow.p - a pickled DataFrame of daily (at noon) flow [m3/s] over many years, e.g. 1980-2022, organized by mean size, with scaling factors applied, which looks like:

```
                        columbia     fraser    squamish  ...  deschutes  nf_skokomish    wilson
1980-01-01 12:00:00          NaN   934.0128         NaN  ...        NaN      6.824360       NaN
1980-01-02 12:00:00          NaN   918.8256         NaN  ...        NaN      6.031488       NaN
1980-01-03 12:00:00          NaN   896.0448         NaN  ...        NaN      5.040399       NaN
1980-01-04 12:00:00          NaN   893.8752         NaN  ...        NaN      4.219210       NaN
1980-01-05 12:00:00          NaN   863.5008         NaN  ...        NaN      3.737824       NaN
...                          ...        ...         ...  ...        ...           ...       ...
2020-12-27 12:00:00  7066.117940  1632.6240  163.058484  ...  16.525168     20.416446  0.433961
2020-12-28 12:00:00  6752.701419  1631.2680  149.448293  ...  15.437212     16.423771  0.398639
2020-12-29 12:00:00  7094.610351  1652.9640  131.500038  ...  14.055214     15.036246  0.351542
2020-12-30 12:00:00  6382.300075  1650.2520  121.563524  ...  15.113766     18.321000  0.501242
2020-12-31 12:00:00  6781.193830  1614.9960  161.358431  ...  32.050593     22.738428  0.649260
```

---

`add_ec_historical.py` an interactive program for adding historical EC rivers, one year at a time, to the DataFrame (++) created by `make_historical.py`.  It exists because make_historical sometimes skips a year due to server timeouts, and it is quicker to fill in important rivers like the Fraser by hand than to run make_historical many times.

---

`plot_historical.py` reads in (++) and makes time series plots for each river (9 per figure).

OUTPUT:

(*)/Data_historical_plots/ALL_flow_plot_[1-6].png

---

`make_historical_temperature.py`  is much like make_historical, but for temperature.

INPUT: like for make_historical, but for EC data we only query the last full year. (Is this because temperature is not available in the older tables?).

OUTPUT:

(++T) = (*)/Data_historical/ALL_temperature.p

**NOTE**: I refactored this to work in pre/river1 on 2023.04.10 and did a bit of testing, but for the output I just copied what we had from pre/river. The reason was that the Fraser temperature from 2022 looked crappy in the summer and I did not want to make it what we rely on in making our climatology (which is all we use for making forcing).

---

`plot_historical_temperature.py` reads in (++T) and makes time series plots.

OUTPUT

- (*)/Data_historical_plots/ALL_temperature_usgs_plot_1.png (currently 12 rivers, 1980-2020, v. gappy)
- (*)/Data_historical_plots/ALL_temperature_ec_plot_1.png (currently 6 rivers, 2020 only, gappy)

---

`make_climatology.py` averages the flow (and temperature if you tell it to) historical records by yearday.  There are no missing values, so these can be used as a backup plan when there is missing data.  Note that only about a third of the rivers have enough temperature data to make climatologies, and that for EC rivers the data is just from the last year.

INPUT:

(++) and (++T)

PLOT OUTPUT:

- (*)/Data_historical_plots/CLIM_flow_plot.png
- (*)/Data_historical_plots/CLIM_temp_plot.png

OUTPUT:

- (*)/Data_historical/CLIM_flow.p
- (*)/Data_historical/CLIM_temp.p

which are pickled Dataframes that look like:

```
        columbia       fraser    squamish     clowhom  ...  dungeness  deschutes  nf_skokomish    wilson
1    7576.033848  1072.712859  166.252107  247.498092  ...  12.801400  19.496954      7.566814  0.434497
2    7497.434094  1076.630927  195.249401  299.487419  ...  13.702930  19.988004      8.041846  0.442161
3    7459.116713  1083.743863  210.613707  303.830145  ...  13.907126  20.403584      8.105456  0.446619
4    7600.596272  1120.439649  209.064770  405.648362  ...  15.416402  18.821637      8.157324  0.451402
5    7586.841314  1087.679571  209.648390  396.917928  ...  16.077416  20.875031      8.282125  0.540191
..           ...          ...         ...         ...  ...        ...        ...           ...       ...
362  6720.410019  1141.341893  182.373620  264.022683  ...  13.188001  22.696259      8.268105  0.499669
363  7120.253521  1130.302068  175.401125  255.510735  ...  12.598818  22.151807      8.219068  0.512987
364  7483.056888  1107.785854  173.989952  215.525272  ...  12.031428  22.184056      7.490773  0.505190
365  7487.805624  1085.838498  163.153582  208.123746  ...  11.523762  22.732303      7.370944  0.467800
366  7525.558068  1205.385382  115.413047  175.890028  ...  13.043897  26.937935     10.197154  0.425753
```

---

`plot_river_map.py` makes a nice plot of the river tracks, and colors the ones with temperature data red.

INPUT:
- (*)/river_info.p - pickled DataFrame of river_info
- (*)/tracks/[river name].p - lon, lat tracks for carving river channels in a grid
- (*)/Data_historical/CLIM_temp.p
- Ldir['grid'] / 'river_info.csv'

PLOT OUTPUT:
- (*)/Data_historical_plots/river_map.png

NOTE: You need to associate a [gridname] with a [ctag] for this to work.
