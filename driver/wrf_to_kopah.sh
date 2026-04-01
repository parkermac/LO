#!/bin/bash

# code to copy today's WRF forecast files to kopah

dstr=`date -u +%Y%m%d`00
indir0=/gscratch/macc/parker/LO_data/wrf/
indir=$indir0$dstr/
# echo $indir
s5cmd sync $indir s3://liveocean-pmacc/LO_data/wrf/$dstr/ > wrf_to_kopah.log &