#!/bin/bash

# Run multiple cast extraction jobs at once

# Example run command:
# ./multi_source_cast_extraction cas6_v0_live

otype=bottle
year=2017

for source in ecology dfo nceiCoastal nceiSalish
do
 python ./extract_casts_fast.py -gtx $1 -source $source -otype $otype -year $year > $source.log &
done