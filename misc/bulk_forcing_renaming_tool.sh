#!/bin/bash

# This custom code will rename a series of forcing folders over
# a date range. For example if you wanted to rename ocnG01 -> ocnG02
# in all folders in LO_output/forcing/cas7 over some date range.

start_date="2024-11-22"
end_date="2024-11-23"

current_date="$start_date"

while [[ "$current_date" <= "$end_date" ]] ; do
  echo "$current_date"
  current_date=$(date --date="$current_date + 1 day" +"%Y-%m-%d")
done