# README for LO/extract/river

#### This code is for extracting time series of river flow and tracer concentration from forcing files from specific runs. The files can then be used in budget calculations, for example in LO/extract/tef2/tracer_budget.py

---

`extract_rivers.py` is the most generic tool. You need to use command line arguments to tell it which [gridname] because this is the primary organizing category for forcing files. Then you also need to tell it which river forcing to use. e.g. [riv] = riv00, and the date range.

It extracts transport and tracer time series, one per day at noon UTC, saved in an xarray Dataset with contents like:

```
<xarray.Dataset>
Dimensions:    (time: 1, riv: 45)
Coordinates:
  * time       (time) datetime64[ns] 2019-07-04T12:00:00
  * riv        (riv) <U13 'clowhom' 'squamish' 'fraser' ... 'umpqua' 'coquille'
Data variables:
    transport  (time, riv) float64 437.2 1.102e+03 4.534e+03 ... 47.75 8.046
    salt       (time, riv) float64 0.0 0.0 0.0 0.0 0.0 ... 0.0 0.0 0.0 0.0 0.0
    temp       (time, riv) float64 8.507 8.507 12.37 16.74 ... 20.15 21.69 21.69
    Oxyg       (time, riv) float64 350.0 350.0 350.0 350.0 ... 350.0 350.0 350.0
    NH4        (time, riv) float64 0.0 0.0 0.0 0.0 0.0 ... 0.0 0.0 0.0 0.0 0.0
    NO3        (time, riv) float64 5.0 5.0 2.737 5.0 5.0 ... 5.0 5.0 5.0 5.0 5.0
    Phyt       (time, riv) float64 0.0 0.0 0.0 0.0 0.0 ... 0.0 0.0 0.0 0.0 0.0
    Zoop       (time, riv) float64 0.0 0.0 0.0 0.0 0.0 ... 0.0 0.0 0.0 0.0 0.0
    SDeN       (time, riv) float64 0.0 0.0 0.0 0.0 0.0 ... 0.0 0.0 0.0 0.0 0.0
    LDeN       (time, riv) float64 0.0 0.0 0.0 0.0 0.0 ... 0.0 0.0 0.0 0.0 0.0
    TIC        (time, riv) float64 300.0 300.0 300.0 300.0 ... 300.0 300.0 300.0
    TAlk       (time, riv) float64 300.0 300.0 300.0 300.0 ... 300.0 300.0 300.0
```
Note: this example just has one time, but in general these are long time series.

The output is saved as NetCDF in:

**LO_output/pre/river1/[gridname]\_[riv]/Data_roms/extraction_[date range].nc**

and you can plot selected results using `plot_rivers.py`.

---

`extract_rivers_cas6_v0_live.py` is just like extract_rivers.py but hard-coded for the full time range of the cas6_v0_live forcing on apogee.

The output is saved as NetCDF in:

**LO_output/pre/river1/lo_base/Data_roms/extraction_2016.12.15_2023.05.21.nc**

and you can plot selected results using `plot_rivers_cas6_v0_live.py`.
