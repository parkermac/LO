#!/bin/bash

# code to push a sequence of files to azure

indir="/dat1/parker/LO_roms/cas7_t0_x4b/averages/"

for (( year=2014; year<=2014; year++)); do
    for (( month=1; month<=2; month++)); do
        if [ $month -lt 10 ]; then
            dstr=$year'_0'$month
        else
            dstr=$year'_'$month
        fi
        fn=$indir'monthly_mean_'$dstr'.nc'
        python ./copy_to_azure.py -fn $fn > monthly_mean.log
    done
done
