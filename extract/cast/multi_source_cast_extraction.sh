#!/bin/bash

# Run multiple cast extraction jobs for many years. Does both bottles and casts, using
# all available sources. Note the use of "_" in naming the log file, which is
# required, otherwise bash thinks I am referencing the variables $source_ and $otype_
# which do not exist.

# The three REQUIRED command line arguments are:
# 1. gtagex [e.g. cas7_t0_x4b]
# 2. roms out number [1-4]
# 3. list_type [average, hourly]

# NOTE there is a problem with using list_type = lowpass because these are missing the
# requested variable "chlorophyll" (and maybe others).

# Currently you have to edit this by hand to change the year(s).

# Example run command:
# ./multi_source_cast_extraction.sh cas7_t0_x4b 2 hourly > x4b_2014.log &

# for year in 2013 2014 2015 2016 2017 2018 2019 2020 2021 2022 2023
for year in 2015
do

    otype=bottle
    # for source in nceiCoastal LineP nceiPNW WOD
    for source in dfo1 ecology_nc nceiSalish nceiCoastal LineP nceiPNW WOD
    do
    python ./extract_casts_fast.py -gtx $1 -ro $2 -lt $3 -source $source -otype $otype -year $year > ./$1"_"$source"_"$otype"_"$year.log
    done

    otype=ctd
    # for source in LineP NHL
    for source in dfo1 ecology_nc LineP NHL ocnms_ctd
    do
    python ./extract_casts_fast.py -gtx $1 -ro $2 -lt $3 -source $source -otype $otype -year $year > ./$1"_"$source"_"$otype"_"$year.log
    done

done