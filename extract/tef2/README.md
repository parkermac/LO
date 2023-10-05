# README for tef2

#### This is new and improved code for doing Total Exchange Flow (TEF) calculations from LO ROMS history files. A few of the changes from the original tef code are that it: (i) works on diagonal and multi-part sections, (ii) automates finding rivers and bounding sections for the segments, and (iii) is built with xarray instead of netCDF4.

The basic idea of TEF is that you extract transport of volume, salt, and other tracers through a bunch of "sections" across an estuarine channel. These are then organized into salinity classes (or bins in a histogram) and the transport in each salinity bin is tidally averaged. So you have a time series of transports in each salinity bin. Then you can further condense the information by adding up all the transport of a given sign and calling these Qin and Qout. Similarly the transport-weighted salinity is called Sin and Sout. The result is that the whole complicated exchange flow through a section is expressed as four terms.

This code also allows you to gather time series of tracer content in the segments between sections, and use these to create tracer budgets for different volumes.

For a complete description of TEF and its calculation please see:
- MacCready, P. (2011). Calculating Estuarine Exchange Flow Using Isohaline Coordinates. Journal of Physical Oceanography, 41, 1116-1124. doi:10.1175/2011JPO4517.1 [Describes the concept and its application to the Columbia River]
- Lorenz, M., Klingbeil, K., MacCready, P., & Burchard, H. (2019). Numerical issues of the Total Exchange Flow (TEF) analysis framework for quantifying estuarine circulation. Ocean Science Discussions, 15, 601-614. doi:10.5194/os-2018-147 [Describes an important improvement to the calculation method, one that is implemented in this code.]

---

#### Workflow Summary

1. Use the GUI tool `create_sections.py` to define the sections.
2. Make by hand `bounding_sections.txt`, a little text file with the names of the open boundary sections, and save it in the same directory that the section definitions went to.
3. Use `plot_collection.py` to look at the results.
4. Run `create_sect_df.py` to turn the section ends you just defined into exact indices on the u- and v-grids.
5. Run `extract_sections.py` for a given model run and date range.
6. Run `process_sections.py` on the output from extract_sections.py.
7. Run `bulk_calc.py` on the output from process_sections.py.
8. Run `bulk_plot.py` to make plots of the results for each section.
9. Run `create_river_info.py` to save point source location information for a run.
10. Run `create_seg_info_dict.py` to generate a dict of information about each of the segments, including all the j,i indices on the rho-grid in each segment.
11. Run `extract_segments.py` on a given model run and date range to get volume-integrals vs. time for tracers in each segment.
12. Run `tracer_budget.py` on a given model run and date range to see if all the extractions you did add up!

#### Example Commands for the Workflow ("run" implies running on mac or pc in ipython)

```
run create_sections -g cas6 -ctag c0
run create_sections -g cas7 -ctag c0

hand edit bounding_sections.txt to have names of the open boundary sections

run plot_collection -gctag cas6_c0
run plot_collection -gctag cas7_c0

run create_sect_df -gctag cas6_c0
run create_sect_df -gctag cas7_c0

python extract_sections.py -gtx cas6_v00_uu0m -ctag c0 -0 2022.01.01 -1 2022.12.31 > extract.log &
[1 hour/year on apogee, salt only]

python extract_sections.py -gtx cas7_trapsV00_meV00 -ro 3 -his_num 1 -ctag c0 -get_bio True -0 2017.01.01 -1 2017.01.10 > trapsV00.log &
[2 min/day on perigee with get_bio True, which would come to 12 hours per year]

[then you could transfer the results to your mac for further processing]

run process_sections.py -gtx cas6_v00_uu0m -ctag c0 -0 2022.01.01 -1 2022.12.31
[10 minutes/year on mac, salt only]

run process_sections.py -gtx cas7_trapsV00_meV00 -ctag c0 -0 2017.01.01 -1 2017.01.10
[6 minutes/10 days on mac with get_bio True]

python bulk_calc.py -gtx cas6_v00_uu0m -ctag c0 -0 2022.01.01 -1 2022.12.31 > bulk.log &
[10 minutes/year on mac, salt only]

python bulk_calc.py -gtx cas7_trapsV00_meV00 -ctag c0 -0 2017.01.01 -1 2017.01.10 > bulk.log &


run create_river_info -gridname cas6 -frc riv00 -dstr 2019.07.04 -test True
run create_river_info -gridname cas7 -frc trapsV00 -dstr 2017.07.04 -test True
[the flags point to a file of river forcing for the run you have copied or created]
[-test True makes a nice plot of point sources]

run create_seg_info_dict -gctag cas6_c0 -riv riv00
run create_seg_info_dict -gctag cas7_c0 -riv trapsV00

python extract_segments.py -gtx cas6_v00_uu0m -ctag c0 -riv riv00 -0 2022.01.01 -1 2022.12.31 > seg_extract.log &
[need to first put LO_output/extract/tef2/seg_info_dict_cas6_c0_riv00.p onto apogee]
python extract_segments.py -gtx cas7_trapsV00_meV00 -ro 3 -his_num 1 -ctag c0 -get_bio True -riv trapsV00 -0 2017.01.01 -1 2017.01.10 > seg_extract.log &
[need to first put LO_output/extract/tef2/seg_info_dict_cas7_c0_trapsV00.p onto apogee]
[1 hour/year on apogee, salt only]

Need to add more about tracer budget...

```
----

## Program Descriptions

---

`create_sections.py` is an interactive graphical tool for defining sections. A single section is saved as a pickled DataFrame with information like:
```
            x          y
0 -122.829951  48.125612
1 -122.822645  48.201052
2 -122.762454  48.228808
```
In this case the section is defined by the longitude (x) and latitude (y) of three points. Two points would be more common. In order to make full use of TEF capabilities such as tracer budgets you should make sure that the sections you define start and end on land.

Sign convention: When we extract transport across a section we define the positive direction in relation to the starting point of the section, using a right hand rule. "o" is the first point on the section line in the sketch below.

```
   +
o-----
   -
```
The output is organized using a gridname and a collection tag [ctag]. For example, my first meaningful collection of sections had [ctag] = c0 and it was on the cas6 grid. This is identified with the combination [gctag] = cas6_c0.

The output ends up in a folder `LO_output/extract/tef2/sections_[gctag]` and each section is a pickled DataFrame named for the section.

**IMPORTANT NOTE**: In order for the segment generation code below to work you need to specify a list of sections that form the open boundaries of your collection. You create this list by hand as a text file. For example, for cas6_c0 I have created `LO_output/extract/tef2/sections_cas6_c0/bounding_sections.txt` which contains:
```
sog7
jdf1
```
---

`plot_collection.py` is a utility program to plot a [gctag] collection of sections so you can look at what you have.

---

`create_sect_df.py` is a tool to turn all the sections whose endpoints are in a [gctag] into a **single** DataFrame called `LO_output/extract/tef2/sect_df_[gctag].p` with contents like:
```
        sn    i     j  irp   jrp  irm   jrm uv  pm
0     ss12  537   652  537   653  537   652  v   1
1     ss12  538   652  538   653  538   652  v   1
2     ss12  539   652  539   653  539   652  v   1
3     ss12  540   652  540   653  540   652  v   1
4     ss12  540   653  540   653  541   653  u  -1
...    ...  ...   ...  ...   ...  ...   ... ..  ..
996    ss7  555   657  555   658  555   657  v   1
997   sog6  134  1227  134  1227  135  1227  u  -1
998   sog6  134  1228  134  1228  135  1228  u  -1
999   ss26  503   658  503   658  504   658  u  -1
1000  ss26  503   658  503   658  503   659  v  -1
```
There is one row for each point. Each point is on either the u- or v-grid. The numbering here is the same as what we get using xarray to read a history file.
- sn tells you which section it is on
- i,j are the indices (column,row) for the the velocity point
- irp,jrp are the indices for the rho-grid point on the positive side
- irm,jrm are the indices for the rho-grid point on the negative side
- uv indicates if the i,j are on the u- or v-grid
- pm is what to multiply the velocity by to get from E/N positive to the section frame of reference.

In principle you should be able to recycle a collection of sections for a different grid when you run this code, but I have not tried.

The creation of the section indices uses a 'stairstep' algorithm that marches along the diagonal direction, finding the u- or v-grid points that stay closest to the desired path. It also trims any points that are on the land mask (mask_rho = 0).

There is some useful graphical output from this program so you should run it on your laptop.

---

`extract_sections.py` does the hard work of extracting transport and tracer values (interpolated from the rho-grid to the u- or v-grid) from a sequence of hourly history files from a ROMS run. Of course the run has to have the same grid that you specified when running `create_sect_df.py`. As usual in the LO system you use command line arguments to tell it which [gtagex], [ctag], time range, and whether or not to get bio variables.

To speed things up this calls `extract_sections_one_time.py` with many simultaneous subprocess calls. This is the code to look at to see which bio variables are being extracted.

The output ends up with the full raw extraction in a NetCDF file, one for each section and named for the section, e.g. 'ai1.nc'. The output, and that of subsequent steps, goes into:

**LO_output/extract/[gtagex]/tef2 = (+)**

in a folder called

**(+)/extractions_[date range]**

Here is an example of how the results for one section are packed in the NetCDF file:
```
<xarray.Dataset>
Dimensions:  (time: 3, p: 8, z: 30)
Coordinates:
  * time     (time) datetime64[ns] 2021-07-04 ... 2021-07-04T02:00:00
Dimensions without coordinates: p, z
Data variables:
    h        (p) float64 17.14 20.32 22.74 21.58 20.65 18.38 16.7 14.99
    dd       (p) float64 ...
    zeta     (time, p) float32 0.2535 0.2531 0.2528 ... 0.0301 0.03066 0.03053
    salt     (time, z, p) float32 ...
    vel      (time, z, p) float64 ...
    DZ       (time, z, p) float64 0.7218 0.864 0.9729 ... 0.2256 0.217 0.2074
```
The dimension "p" means a point on the stairstep section. "dd" is the point width [m], and "DZ" [m] is the vertical thickness of each cell.

Except when you are testing you have to run this code on the remote computer where the ROMS history files are stored. Since you created your sect_df_[gctag].p on your laptop then you have to copy it by hand to the remote machine.

---

`process_sections.py` goes through all the extracted sections and does the salinity binning, saving the output as NetCDF xarray Datasets (one for each section) in the folder

**(+)/processed_[date range]**

This is fast enough to run on your laptop, after you have copied the results of `extract_sections.py`.

Warning: if you look around line 90 you will see that this is hard-coded to use 1000 salinity bins between 0 and 36.

It automatically figures out which data variables were extracted, and also processes salt-squared, for variance budgets.

Here is an example of what is in the Datasets:
```
<xarray.Dataset>
Dimensions:  (time: 8761, sbins: 36)
Coordinates:
  * time     (time) datetime64[ns] 2022-01-01 2022-01-01T01:00:00 ... 2023-01-01
  * sbins    (sbins) float64 0.5 1.5 2.5 3.5 4.5 ... 31.5 32.5 33.5 34.5 35.5
Data variables:
    qnet     (time) float64 1.291e+06 2.207e+06 ... 2.059e+06 2.289e+06
    fnet     (time) float64 -3.724e+09 -1.807e+10 ... -5.21e+09 -1.007e+10
    ssh      (time) float64 -0.2873 -0.8149 -1.287 ... 0.05648 -0.2519 -0.4378
    salt     (time, sbins) float64 0.0 0.0 0.0 0.0 0.0 ... 1.811e+07 0.0 0.0 0.0
    q        (time, sbins) float64 0.0 0.0 0.0 0.0 0.0 ... 5.62e+05 0.0 0.0 0.0
    salt2    (time, sbins) float64 0.0 0.0 0.0 0.0 0.0 ... 5.834e+08 0.0 0.0 0.0
```
This is from running with -test True, which decreases the number os salinity bins to 36, whereas in production runs we use 1000.

---

`bulk_calc.py` boils down the results of process_sections.py into some number of 'in' and 'out' layers, using code from Marvin Lorenz. Each of these is a time series, tidally averaged and subsampled to daily at noon UTC of each day. The results end up in xarray Datasets, saved as NetCDF, one for each section in the collection. As a convenience we now include 'qprism' as a time series [m3 s-1].

**(+)/bulk_[date range]**

There can be more than two layers!

Here is what is in the Datasets:
```
<xarray.Dataset>
Dimensions:  (time: 363, layer: 30)
Coordinates:
  * time     (time) datetime64[ns] 2022-01-02T12:00:00 ... 2022-12-30T12:00:00
  * layer    (layer) int64 0 1 2 3 4 5 6 7 8 9 ... 20 21 22 23 24 25 26 27 28 29
Data variables:
    salt     (time, layer) float64 31.22 33.0 nan nan nan ... nan nan nan nan
    q        (time, layer) float64 9.558e+04 -1.171e+05 nan nan ... nan nan nan
    salt2    (time, layer) float64 974.8 1.089e+03 nan nan ... nan nan nan nan
    qnet     (time) float64 -2.153e+04 -1.27e+04 ... 4.075e+03 5.496e+03
    fnet     (time) float64 -1.019e+10 -1.064e+10 ... -4.705e+09 -3.847e+09
    ssh      (time) float64 0.2146 0.2581 0.2727 0.2139 ... 0.2784 0.2749 0.2785
    qabs     (time) float64 7.63e+05 1.493e+06 1.536e+06 ... 1.151e+06 5.106e+05
    qprism   (time) float64 3.815e+05 7.466e+05 ... 5.756e+05 2.553e+05
```

---

`bulk_plot.py` is for making plots of all the bulk fields, showing both the multi-layer results of `bulk_calc.py` and the results of pushing them into two layers.

The results are a bunch of png's, one for each section in

**(+)/bulk_plots_[date range]**

---

#### Next there are a series of steps to define and extract data for _segments_, which are the volumes between sections. These are essential to form the storage term in volume-integrated budgets.

---

`create_river_info.py` processes a user-specified river forcing file to create and save a DataFrame that has info about each river. Here is what it looks like for a cas6 file for [riv] = traps2:
```
               name     ii      jj  dir  sgn     iu      ju     iv     jv   irho    jrho
0           clowhom  453.0  1169.0  0.0 -1.0  452.0  1169.0    NaN    NaN  452.0  1169.0
1          squamish  484.0  1169.0  0.0  1.0  483.0  1169.0    NaN    NaN  484.0  1169.0
2            fraser  641.0  1082.0  0.0 -1.0  640.0  1082.0    NaN    NaN  640.0  1082.0
3            tsolum  250.0  1160.0  0.0 -1.0  249.0  1160.0    NaN    NaN  249.0  1160.0
4            oyster  227.0  1179.0  0.0  1.0  226.0  1179.0    NaN    NaN  227.0  1179.0
..              ...    ...     ...  ...  ...    ...     ...    ...    ...    ...     ...
268       Hartstene  536.0   675.0  0.0  1.0  535.0   675.0    NaN    NaN  536.0   675.0
269  Seashore Villa  529.0   631.0  1.0  1.0    NaN     NaN  529.0  630.0  529.0   631.0
270            LOTT  528.0   622.0  1.0  1.0    NaN     NaN  528.0  621.0  528.0   622.0
271      Rustlewood  525.0   671.0  0.0  1.0  524.0   671.0    NaN    NaN  525.0   671.0
272         Shelton  503.0   655.0  1.0  1.0    NaN     NaN  503.0  654.0  503.0   655.0
```
The main job of this code is to generate a list of indices on the rho grid (jrho, irho) that each river flows INTO. It should work for sources on the u and v grids, and also for vertical sources on the rho grid, although I have not tested that yet.

The rho-grid indices are then used in `create_seg_info_dict.py` when it is assembling a list of rivers in each segment.

If you run the code with -test True it makes a useful map plot which you can use to see visually that all the sources look good.

The DataFrame is pickled and saved as (for the example above):

**LO_output/extract/tef2/riv_df_[gridname]_[riv].p**

e.g. riv_df_cas6_traps2.p

---

`create_seg_info_dict.py` is some rather involved (but useful) code to define tef2 segments. Mainly we get all the j,i indices on the rho grid for all the segments between sections. For each segment we also get a list of the bounding sections and the rivers.

The results are packed in a dict of dicts called seg_info_dict, and saved as a pickled dict in:

**LO_output/extract/tef2/seg_info_dict_[gctag]_[riv].p**

e.g. seg_info_dict_cas6_c0_traps2.p

seg_info_dict has one key for each segment, e.g. 'mb8_m' which is a section name and a sign (plus or minus side - see the section sign convention above). Note that it is sort of random which section gets assigned to be the key defining the segment - it could be any of the bounding sections. One advantage of this system is that it is all done automatically, and so we don't have to make a custom list like in LO/tef.

Then seg_info_dict['mb8_m'] is a dict with three keys:
['ji_list', 'sns_list', 'riv_list']
- ji_list is a list of tuples (j,i) that are indices of points in the rho grid in the segment
- sns_list is a list of all the bounding sections+signs (sns is short for section name + sign), including their signs. Note that the segment key is one member of this list
- riv_list is a list of the rivers or point sources entering the segment

It also creates and saves a pickled pandas DataFrame with columns:
['volume m3', 'area m2', 'lon', 'lat']
which are the volume (with zeta = 0), surface area, and mean lon and lat of each segment.

**LO_output/extract/tef2/vol_df_[gctag].p**

---

`extract_segments.py` extracts integrated tracer values as hourly time series for each of the segments. The results are saved in an xarray Dataset:

**(+)/segments_[date range]\_[gctag]\_[riv].nc**

Uses `extract_segments_one_time.py` and subprocess to speed execution.

---

`tracer_budget.py` forms volume-integrated budgets, combining section, segment, and river extractions.

NOTE: To use this for a case with rivers you need to have extracted river time series from the appropriate forcing files, using `LO/extract/river/extract_rivers.py`. Also the time range of the river extraction has to match the time range of the section and segment extractions.

UNDER CONSTRUCTION

---
