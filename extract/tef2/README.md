# README for tef2

#### This is new and improved code for doing Total Exchange Flow (TEF) calculations from LO ROMS history files. A few of the changes from the original tef code are that it: (i) works on diagonal and multi-part sections, (ii) automates finding rivers and bounding sections for the segments, and (iii) is built with xarray instead of netCDF4.

The basic idea of TEF is that you extract transport of volume, salt, and other tracers through a bunch of "sections" across an estuarine channel. These are then organized into salinity classes (or bins in a histogram) and the transport in each salinity bin is tidally averaged. So you have a time series of transports in each salinity bin. Then you can further condense the information by adding up all the transport of a given sign and calling these Qin and Qout. Similarly the transport-weighted salinity is called Sin and Sout. The result is that the whole complicated exchange flow through a section is expressed as four terms.

For a complete description please see: MacCready, P. (2011). Calculating Estuarine Exchange Flow Using Isohaline Coordinates. Journal of Physical Oceanography, 41, 1116-1124. doi:10.1175/2011JPO4517.1

---

#### Workflow Summary

[need to add this]

----

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
The output is organized using a gridname and a collection tag [ctag]. For example, my first meaningful collection of sections had [ctag] = c0 and it was on the cas6 grid. This is identified with combination [gctag] = cas6_c0.

The output ends up in a folder `LO_output/extract/tef2/sections_[gctag]` and each section is a pickled DataFrame named for the section.

---

`plot_collection.py` is a utility program to plot a [cgtag] collection of sections so you can look at what you have.

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
There is one row for each point. Each point is on either the u- or v-grid. The numbering here is the same in what we get using xarray to read a history file.
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

To speed things up this calls `get_one_section.py` with many simultaneous subprocess calls. This is the code to look at to see which bio variables are being extracted.

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

`process_sections.py` goes through all the extracted sections and does the salinity binning, saving the output in picked dicts of numpy arrays in the folder

**(+)/processed_[date range]**

This is fast enough to run on your laptop, after you have copied the results of `extract_sections.py`.

Warning: if you look around line 90 you will see that this is hard-coded to use 1000 salinity bins between 0 and 36.

By now you may have noticed that I jump around aimlessly, sometimes saving output in picked DataFrames, picked dicts, or NetCDF. Sorry. Often the choice is for convenience, or to smooth reuse of old code from LO/tef. It also is a result of the incredible speed advantages offered by nco tools like ncrcat.

---

`bulk_calc.py` boils down the results of process_sections.py into some number of 'in' and 'out' layers, using code from Marvin Lorenz. Each of these is a time series, tidally averaged and subsampled to daily at noon UTC of each day. The results end up in pickled dicts, one for each section in

**(+)/bulk_[date range]**

There can be more than two layers!

---

`bulk_plot.py` is for making plots of all the bulk fields, showing both the multi-layer results of `bulk_calc.py` and the results of pushing them into two layers.

The results are a bunch of png's, one for each section in

**(+)/bulk_plots_[date range]**
