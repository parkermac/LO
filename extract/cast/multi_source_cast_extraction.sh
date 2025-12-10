#!/bin/bash

# Run multiple cast extraction jobs for one or more years. Does both bottles and casts, using
# all available sources. Note the use of "_" in naming the log file, which is
# required, otherwise bash thinks I am referencing the variables $source_ and $otype_
# which do not exist.

# The five REQUIRED command line arguments are:
# 1. gtagex [e.g. cas7_t0_x4b]
# 2. roms out number [1-4]
# 3. list_type [average, hourly]
# 4. starting year [e.g. 2017]
# 5. ending year [e.g. 2020; use 2017 to just get one year]

# NOTE there is a problem with using list_type = lowpass because these are missing the
# requested variable "chlorophyll" (and maybe others).


# Example run command:
# ./multi_source_cast_extraction.sh cas7_t0_x4b 2 hourly 2016 2016 > x4b_2016.log &
# ./multi_source_cast_extraction.sh cas7_t1_x11ab 1 average 2016 2016 > x11ab_2016.log &

gtx=$1
ro=$2
lt=$3
year0=$4
year1=$5

for (( year = $year0; year <= $year1; year++ ))
do

    otype=bottle
    for source in dfo1 ecology_nc nceiSalish nceiCoastal LineP nceiPNW WOD
    do
    echo "Executed ./extract_casts_fast.py -gtx $gtx -ro $ro -lt $lt -source $source -otype $otype -year $year > ./$gtx"_"$source"_"$otype"_"$year.log"
    python ./extract_casts_fast.py -gtx $gtx -ro $ro -lt $lt -source $source -otype $otype -year $year > ./$gtx"_"$source"_"$otype"_"$year.log
    done

    otype=ctd
    for source in dfo1 ecology_nc LineP NHL ocnms_ctd
    do
    echo "Executed ./extract_casts_fast.py -gtx $gtx -ro $ro -lt $lt -source $source -otype $otype -year $year > ./$gtx"_"$source"_"$otype"_"$year.log"
    python ./extract_casts_fast.py -gtx $gtx -ro $ro -lt $lt -source $source -otype $otype -year $year > ./$gtx"_"$source"_"$otype"_"$year.log
    done

done