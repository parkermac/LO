Notes from Charles Seaton for the CRITFC extraction.
2020.11.10

Here is the python code and SELFE files needed to generate nudging files for the CRITFC-CMOP forecast. 

1) makenudge.sh is a bash shell script that sets up the python environment variables and the python script arguments, which will need to be modified for your environment.

2) gen_cmop_nudge.py is the python script, it is written for python 2.7 and uses the python libraries:
numpy, netCDF4 , scipy
It takes 6 arguments (examples shown in the makenudge.sh file):
the full path to the SELFE hgrid.ll file
the full path to the SELFE vgrid.in file
the full path to the LiveOcean ocean_depth file containing the variables: rho_lon, rho_lat, and depths
the base path to the LiveOcean daily output directory
the full path to the directory where output files will be written
the run date in YYYY-MM-DD format

3/4) hgrid.ll and vgrid.in are the SELFE horizontal and vertical grid files required to generate the nudging file.

Please let me know if anything is unclear or if you have any problems working with these files.

Thanks again for providing us with your LiveOcean output for our forecast nudging and for offering to set up generation of the nudging files on your systems.

Charles