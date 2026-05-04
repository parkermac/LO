#!/bin/bash

# This script allows you to send a file to kopah for sharing
# using s5cmd.

# Complete path or relative path of the input file
in_fn=$1

s5cmd cp --acl public-read $1 s3://pm-share/
echo "URL = https://s3.kopah.uw.edu/pm-share/[file name]"