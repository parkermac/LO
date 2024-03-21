#### OCNMS mooring data in LO format:
Quick description: variables for the 10 current, ocnms moorings are processed in LO format for years 2011 - 2023. 
these represent hourly averages of 10 minute mooring data. Saved variables include: SA, SP, IT, CT, DO, SIG0, and P. 

#### OCNMS nomenclature + site designations: 
Moored sensors were first deployed in 2000 and 10** lightweight, coastal moorings are currently 
maintained from Cape Elizabeth (“CE”) to Makah Bay (“MB”) at depths of ~15m, ~42m, or ~27m.

North to south sites:
MB = Makah Bay 
CA = Cape Alava
TH = Teahwhit Head
KL = Kalaloch
CE = Cape Elizabeth

The mooring depth is listed with its site id, e.g., CE042 is placed at ~42m; 
KL027 at ~27m; and CE015 is placed at ~15m.

##### OCNMS job_names:
'OCNMS_moorings_current':
**A list of the 10 moorings currently maintained in the sanctuary:
	MB015
	MB042
	CA015
	CA042
	TH015
	TH042
	KL015
	KL027
	CE015
	CE042

'OCNMS_moorings_historical':
A full list of OCNMS moorings. 
Historically, moorings were also deployed at depths >42m. 
These historical mooorings* were not deployed during LO model years 2017 - (2023). 
Or during the long hindcast (2012 - present). 
	MB015
	MB042
	CA015
	CA042
	CA065*
	CA100*
	TH015
	TH042
	TH065*
	KL015
	KL027
	KL050*
	CE015
	CE042
	CE065* 
    
##### OCNMS saved variables:
Part 1: Raw to hourly data: 
Matlab processing steps (OCNMS data are provided as matfiles)
- Raw data saved here are from matfiles saved/shared on a google drive provided by Brandy Cervantes/OSU. 
latest version + access date was 8 December 2023 
- From this source we removed sites that predate 2011 
- Raw data were used to calculate hourly averages for SP, IT, DO and P 
- ~10 days of mooring data were flagged suspect at the end of a 2021 deployment at MB042 and are not included. 
- Several DO values when reading ~slight negative (e.g. -0.000X mg/L) were adjusted to 0 mg/L. (we might want to save a flag variable with the LO files)
- Depths were estimated based on altitude. See OCNMS mooring descriptions online at: 
https://olympiccoast.noaa.gov/science/oceanographic-moorings/data.html
- hourly matfiles are saved here: .../LO_data/obs/ocnms/ocnms_mooring/datasets 
with units: 
time in UTC 
oxy mg/l 
temp deg C insitu
sal practical 
Note: BC relayed that oxy = mgL and time UTC 

Part 2: Processing to LO format: 
- process_moorings.py will process hourly matfiles; calculate SA and SIG0 and convert DO to uM; and 
then save a netcdf with LO formatted data here: ... /LO_output/obs/ocnms/moor
- this provides variables for the 10 current moorings in LO format for years 2011 - 2023. 
- saved variables include hourly data for: SA, SP, IT, CT, DO, SIG0, and P. 

- plot_processed_locations.py will create a figure illustrating data availability for the 10 current moorings


 
