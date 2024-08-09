#!/bin/bash

year=$1
ymd0=$year".07.26"
ymd1=$year".07.31"
gridname=cas7
python ./driver_forcing3.py -g $gridname -0 $ymd0 -1 $ymd1 -f ocn01 > ocn01_$year.log &
python ./driver_forcing3.py -g $gridname -0 $ymd0 -1 $ymd1 -f atm00 > atm00_$year.log &
python ./driver_forcing3.py -g $gridname -0 $ymd0 -1 $ymd1 -tP trapsP00 -f trapsF00 > trapsF00_$year.log &
python ./driver_forcing3.py -g $gridname -0 $ymd0 -1 $ymd1 -f tide00 > tide00_$year.log &
