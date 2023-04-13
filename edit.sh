#!/bin/bash

counter=1

while true; do
    num=$((counter+848258))
    sed -i "12s/\(.\{16\}\).*/\1 848258 $num/" run1.mac
    ((counter++))
    sleep 1
done
