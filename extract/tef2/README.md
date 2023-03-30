# README for tef2

#### This is new and improved code for doing Total Exchange Flow (TEF) calculations from LO ROMS history files. A few of the changes from the original tef code are that it: (i) works on diagonal and multi-part sections, (ii) automates finding rivers and bounding sections for the segments, and (iii) is built with xarray instead of netCDF4.

The basic idea of TEF is that you extract transport of volume, salt, and other tracers through a bunch of "sections" across an estuarine channel. These are then organized into salinity classes (or bins in a histogram) and the transport in each salinity bin is tidally averaged. So you have a time series of transports in each salinity bin. Then you can further condense the information by adding up all the transport of a given sign and calling these Qin and Qout. Similarly the transport-weighted salinity is called Sin and Sout. The result is that the whole complicated exchange flow through a section is expressed as four terms.

For a complete description please see: MacCready, P. (2011). Calculating Estuarine Exchange Flow Using Isohaline Coordinates. Journal of Physical Oceanography, 41, 1116-1124. doi:10.1175/2011JPO4517.1

---

#### The code...

[need a summary of the workflow]

----

`create_sections.py` is an interactive graphical tool for defining sections. A single section is saved as a pickled DataFrame with information like:
```
            x          y
0 -122.829951  48.125612
1 -122.822645  48.201052
2 -122.762454  48.228808
```
In this case the section is defined by the longitude (x) and latitude (y) of three points, although two points would be more common. In order to make full use of TEF capabilities such as tracer budgets you should make sure that the sections you define start and end on land.

Sign convention: When we extract transport across a section we define the positive direction in relation to the starting point of the section, using a right hand rule. "o" is the first point on the section line in the sketch below.

```
   +
o-----
   -
```
