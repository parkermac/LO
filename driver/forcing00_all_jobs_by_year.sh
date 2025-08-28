#!/bin/bash

year=$1
ymd0=$year".07.26"
ymd1=$year".07.31"
gridname=cas7
python ./driver_forcing00.py -g $gridname -0 $ymd0 -1 $ymd1 -do_bio True -f ocnG01 > ocn01_$year.log &
python ./driver_forcing00.py -g $gridname -0 $ymd0 -1 $ymd1 -f atm01 > atm01_$year.log &
python ./driver_forcing00.py -g $gridname -0 $ymd0 -1 $ymd1 -tP trapsP01 -f trapsN00 > trapsN00_$year.log &
python ./driver_forcing00.py -g $gridname -0 $ymd0 -1 $ymd1 -f tide01 > tide01_$year.log &
