#!/bin/bash

counter=1

while true
do
  # execute command here
  ./SLACtut run1.mac
  num=$((counter+848258))
  sed -i "12s/\(.\{16\}\).*/\1 848258 $num/" run1.mac

  # rename files and move them to ../data_raw directory
  mv Tutorial_nt_Tutorial.csv ../data_raw/ntuple_$counter.csv
  mv Tutorial_h2_CHC1_XY.csv ../data_raw/histo_$counter.csv
  ((counter++))

done
