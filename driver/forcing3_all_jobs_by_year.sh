#!/bin/bash

year=$1
ymd0=$year.01.01
ymd1=$year.12.31
gridname=cas7
python ./driver_forcing3.py -g $gridname -0 $ymd0 -1 $ymd1 -f ocn01 > ocn01_$year.log &
python ./driver_forcing3.py -g $gridname -0 $ymd0 -1 $ymd1 -f atm00 > atm00_$year.log &
python ./driver_forcing3.py -g $gridname -0 $ymd0 -1 $ymd1 -f traps00 > traps00_$year.log &
python ./driver_forcing3.py -g $gridname -0 $ymd0 -1 $ymd1 -f tide00 > tide00_$year.log &
