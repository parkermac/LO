#!/bin/bash

# Run multiple cast extraction jobs at once

# Example run command:
# ./multi_source_cast_extraction.sh cas6_v0_live 0 [my run on perigee]
# ./multi_source_cast_extraction.sh cas7_trapsV00_meV00 3 [Aurora's run on perigee]

otype=bottle
year=2017
for source in ecology dfo1 nceiCoastal nceiSalish
# for source in dfo1
do
 python ./extract_casts_fast.py -gtx $1 -ro $2 -source $source -otype $otype -year $year > $source_$otype_$year.log &
done

otype=ctd
year=2017
for source in dfo1 ecology
do
 python ./extract_casts_fast.py -gtx $1 -ro $2 -source $source -otype $otype -year $year > $source_$otype_$year.log &
done