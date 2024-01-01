#!/bin/bash

# Run multiple cast extraction jobs at once, does both bottles and casts, using
# all available sources. Note the use of "_" in naming the log file, which is
# required, otherwise bash thinks I am referencing the variables $source_ and $otype_
# which do not exist.

# Example run commands:
# ./multi_source_cast_extraction.sh cas6_v0_live 0 2017 [my run on perigee]
# ./multi_source_cast_extraction.sh cas7_trapsV00_meV00 3 2017 [Aurora's run on perigee]
# ./multi_source_cast_extraction.sh cas7_t0_x4b 0 2013 [Long Hindcast on apogee]

# The three required command line arguments are:
# 1. gtagex
# 2. roms out number
# 3. year

otype=bottle
for source in ecology dfo1 nceiCoastal nceiSalish
do
 python ./extract_casts_fast.py -gtx $1 -ro $2 -source $source -otype $otype -year $3 > ./$1"_"$source"_"$otype"_"$3.log &
done

otype=ctd
for source in dfo1 ecology
do
 python ./extract_casts_fast.py -gtx $1 -ro $2 -source $source -otype $otype -year $3 > ./$1"_"$source"_"$otype"_"$3.log &
done