#!/bin/bash

year=2017
python ./driver_forcing3.py -g cas6 -0 $year.01.01 -1 $year.12.31 -f ocn00 > ocn00_$year.log &
python ./driver_forcing3.py -g cas6 -0 $year.01.01 -1 $year.12.31 -f atm00 > atm00_$year.log &
python ./driver_forcing3.py -g cas6 -0 $year.01.01 -1 $year.12.31 -f riv00 > riv00_$year.log &
python ./driver_forcing3.py -g cas6 -0 $year.01.01 -1 $year.12.31 -f tide00 > tide00_$year.log &
