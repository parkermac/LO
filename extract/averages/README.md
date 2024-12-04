# README for LO/extract/averages

## Code for creating and plotting monthly averages and climatologies of ROMS output

### **extract_month_mean.py** works on the daily lowpassed.nc files that must have already been created. It averages them into monthly means, with the same data_vars.

- The output files are put in LO_roms/[gtagex]/averages/monthly_mean_[YYYY]_[MM].nc

- **month_mean_worker.py** is a subprocess called by extract_month_mean that does part of the averaging. It allows us to speed up the averaging by working in parallel.

### **climatology_by_month.py** works on the output of extract_month_mean, for example averaging all of the Januarys for all available years into a single "climatological" January.

- The output files are put in LO_roms/[gtagex]/climatologies/monthly_clim_[MM].nc

### **month_clim_grid_plot.py** is a tool to create single plots that summarize the results of making monthly means and climatologies.

- you use command line arguments to specity which variable to use (surface temp or bottom oxygen) and whether to plot the monthly mean or the monthly anomaly (monthly mean minus climatology)

- output is saved as a png in LO_output/plots

- NOTE: you can also make plots of individual months or movies using pan_plot.py with -lt monthly_mean and -pt P_monthly_mean (this plot_type has to be hand-edited to work on temp or oxygen).