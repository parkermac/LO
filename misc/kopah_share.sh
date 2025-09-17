#!/bin/bash

# This script allows you to send a file to kopah for sharing.

# Testing:
# echo "hi" > test.txt
# ./kopah_share.sh test.txt test.txt
# or for a bigger file:
# ./kopah_share.sh test.txt test.txt > test.log &

# Complete path or relative path of the input file
in_fn=$1

# Name of the output file
out_fn=$2

s3cmd put --acl-public $1 s3://pm-share/$2
echo "URL = https://s3.kopah.uw.edu/pm-share/"$2