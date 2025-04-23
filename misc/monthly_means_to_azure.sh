#!/bin/bash

# code to push a sequence of files to azure
# takes about a minute per file and there are 120 files

indir="/dat1/parker/LO_roms/cas7_t0_x4b/averages/"

for (( year=2014; year<=2023; year++)); do
    for (( month=1; month<=12; month++)); do
        if [ $month -lt 10 ]; then
            dstr=$year'_0'$month
        else
            dstr=$year'_'$month
        fi
        fn=$indir'monthly_mean_'$dstr'.nc'
        echo $fn
        python ./copy_to_azure.py -fn $fn > monthly_mean.log
    done
done
