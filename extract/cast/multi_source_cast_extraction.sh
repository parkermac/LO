#!/bin/bash

# Run multiple cast extraction jobs at once

# Example run command:
# ./multi_source_cast_extraction.sh cas6_v0_live

otype=bottle
year=2017
for source in ecology dfo1 nceiCoastal nceiSalish
# for source in dfo1
do
 python ./extract_casts_fast.py -gtx $1 -source $source -otype $otype -year $year > $source_$otype_$year.log &
done

otype=ctd
year=2017
for source in dfo1 ecology
do
 python ./extract_casts_fast.py -gtx $1 -source $source -otype $otype -year $year > $source_$otype_$year.log &
done