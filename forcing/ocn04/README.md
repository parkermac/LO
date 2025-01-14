# README for LO/forcing/ocn04

## Notes from Michael McDonald in an email of 12/30/2024. These are the basis of the edits we will make to turn ocn03 into ocn04.

Parker,
Ah! I think that I see your issue now. The time index will change as new data is added and old data is removed. This is the case for the FMRC "best" datasets, an aggregation of multiple runs and the base "hours since" will change day to day. So, do not use (guess) index values to determine the "date/time" you are seeking, let ncks do this for you. You can use indices for the lat lon (if this works for you), but you can also use the nicer lat/lon decimal equivalent as ncks also converts/finds this in the lat/lon array automatically. 

i.e., *Do not use index values in your ncks queries for time*.

e.g.
NO
-d time,44,100

YES (format is "YYYY-MM-DDThh:mm:ssZ"), -d time,min[,max][,stride]
-d time,2024-12-30T12:00:00Z,2024-12-31T12:00:00Z

e.g. to get 24 hours from  2024-12-30T12:00:00Z to 2024-12-31T12:00:00Z

```
ncks -D5 -v salinity -d time,2024-12-30T12:00:00Z,2024-12-31T12:00:00Z -d lat,39.0,53.0 -d lon,228.95,238.96 "http://tds.hycom.org/thredds/dodsC/FMRC_ESPC-D-V02_s3z/FMRC_ESPC-D-V02_s3z_best.ncd" -O parker_time_test2.nc

```
You will see in the ncks output (with some debugging enabled -D5) like this,

ncks: INFO nco_cln_clc_dff()() reports conversion between systems "s@2024-12-30T00:00:00Z" and "hours since 2024-12-21 12:00:00.000 UTC" is 204.000000
ncks: INFO nco_cln_clc_dff()() reports conversion between systems "s@2024-12-30T23:00:00Z" and "hours since 2024-12-21 12:00:00.000 UTC" is 227.000000

ncks: nco_cln_clc_tm() reports unt_sng=2024-12-30T23:00:00Z bs_sng=hours since 2024-12-21 12:00:00.000 UTC

