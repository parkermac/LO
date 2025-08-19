# README for ecology_nc

These files are WA Dept. of Ecology water column samples, both bottle and CTD. This NetCDF format is a new development as of May 2024. It replaces the "ecology" source.

These were downloaded from:
https://ecology.wa.gov/research-data/monitoring-assessment/puget-sound-and-marine-monitoring/water-column-data
on 2024.05.09. Look for the drop down menu for "Yearly marine water column profiles and nutrient data in netCDF for large data users." The units and other metadata are in the NetCDF files, and on the webpage.

Files go from 1999-2023 (and presumably a new one will be added each year). No carbon data. The processing of other fields was much easier than in the past using Excel files from Ecology. I assumed that "Fluorescence (mg/m^3)" was the same as what we call "Chl (mg m-3)" in the LO format. I also converted DO from mg/L to uM.

The NetCDF files have all the CTD and bottle data in one listing - at depths where there are bottle samples that row has more data. You can explore the data structure by reading, and running, the process_test.py code here. The main thing to be aware of is that there are three dimensions: stations, profiles, and obs, and numerical keys that allow you to associate a given row of obs (one depth of a cast) with a unique profile number, and then associate that profile number with a unique station. This system eliminates the guesswork needed in the past to ensure that we know what a cast was.

The code process_data.py creates both ctd and bottle processed files from a single input file. It assumes a row is from a bottle if the NO3 field has data.

UPDATE: on 2025.08.19 I downloaded and processed the 2024 data. All appeared to go smoothly. Suzan Pool also updated all the NetCDF files but said that the only important edits were for 2005-2007, which I also downloaded and re-processed. She said: "Corrected units for DO concentrations and saturation at SJF000, SJF001, and SJF002 from Sep 2005 to Jun 2007.  They were previously in ml/l, I corrected them to mg/l, and Christopher and I reviewed to confirm success and validity."