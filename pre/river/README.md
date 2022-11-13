# README for the pre/river code

### This code is for gathering and processing river data, and one of the uses for the data is forcing ROMS hindcasts.

### It relies on functions in the module `LO/lo_tools/lo_tools/river_functions.py`, which go out and get data from external sources: USGS, NWS Forecasts, and Environment Canada websites.  The reason for putting these in lo_tools is that they are also used by the forecast system.

#### _NOTE: This was rewritten in November 2022 to make it easier to add one new year to the historical transport data. You have to edit things in make_historical.py each time you do this, and it has to match up with what you look for when making the forcing files, like with riv00._

---

#### CURRENT OUTPUT STATE, 2022.11.13
cas6_v3 has:
- historical flow and flow climatology for 1980-2021
- historical temperature and temperature climatology for 1980-2020

---

The output goes into:

(*) = LO_output/pre/river/[gtag]

where [gtag] could be something like "cas6_v3".  There is no perfect way to name things, but a given collection of rivers is often associated with a specific grid and forcing tag.

2022.10.18 For river forcing code like riv00 we have this hardwired, so that is a bit non-general. Also we have moved away from the [gtag] as an organizing concept. In the future this output directory name could be something more meaningful, like a collection name and a date.

---

`make_river_info.py` organizes river info from different sources.

INPUT:

Information from Sarah Gidding, Ecology, and other places.

OUTPUT:

- (*)/river_info.csv - names, gauges, scaling factors
- (*)/tracks/[river name].p - lon, lat tracks for carving river channels in a grid

---

`make_historical.py` - the main tool for getting historical flow data. This was recoded in November 2022 to make it easier to just do one year at a time and then concatenate it with previous saved data. Much faster and more reliable, but requires you to edit the code each time, and to get up-to-date ROMS extractions first.

INPUT:

Raw data from XML (USGS), or scraped from html (EC), using the functions in river_functions.py.  For the USGS data you can simply hand it a date range and it will quickly return 40 years of daily flow data.  For EC data it is more complicated.  Typically over a range in the last 18 months you can get data using a date range, just like for USGS.  Earlier than this you need to scrape the "historical" data from a table on a web page, a year at a time, but these only exist up to two years ago, e.g. 2020 when it is 2022. If yo do this after mid-year there is a gap in the EC data, e.g. doing this in November 2022 there is a gap for the early months of 2021. A workaround is to use as-run river extractions from LiveOcean forcing files, which are created using code in LO/extract/river.

OUTPUT:

(++) = (*)/Data_historical/ALL_flow_[year0]_[year1].p - a pickled DataFrame of daily (at noon) flow [m3/s] over many years, e.g. 1980-2020, organized by mean size, with scaling factors applied, which looks like:

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

`get_all_flow_df.py` a convenience function for loading the DataFrame (++) created by `make_historical.py`

---

`add_ec_historical.py` an interactive program for adding historical EC rivers, one year at a time, to the DataFrame (++) created by `make_historical.py`.  It exists because make_historical sometimes skips a year due to server timeouts, and it is quicker to fill in important rivers like the Fraser by hand than to run make_historical many times.

---

`plot_historical.py` reads in (++) and makes time series plots for each river (9 per figure).

OUTPUT:

(*)/Data_historical/ALL_flow_plot_[1-6].png

---

`make_historical_temperature.py`  is much like make_historical, but for temperature.

INPUT: like for make_historical, but for EC data we only query the last year.

OUTPUT:

(++T) = (*)/Data_historical/ALL_temperature_[year0]_[year1].p

---

`plot_historical_temperature.py` reads in (++T) and makes time series plots.

OUTPUT

- (*)/Data_historical/ALL_temperature_usgs_plot_1.png (currently 12 rivers)
- (*)/Data_historical/ALL_temperature_ec_plot_1.png (currently 6 rivers)

---

`make_climatology.py` averages the flow (and temperature if you tell it to) historical records by yearday.  There are no missing values, so these can be used as a backup plan when there is missing data.  Note that only about a third of the rivers have enough temperature data to make climatologies, and that for EC rivers the data is just from the last year.

INPUT:

(++) and (++T)

PLOT OUTPUT:

- (*)/Data_historical/CLIM_flow_plot.png
- (*)/Data_historical/CLIM_temp_plot.png

OUTPUT:

- (*)/Data_historical/CLIM_flow_[year0]_[year1].p
- (*)/Data_historical/CLIM_temp_[year0]_[year1].p

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
