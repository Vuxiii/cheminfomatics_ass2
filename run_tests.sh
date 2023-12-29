#!/bin/bash

for i in {0..19}; do
    output=$(mod -f test.py -e "options[$i].func()" 2> /dev/null)
    
    n=$(echo "$output" | grep -oP 'n: \K\d+')
    t=$(echo "$output" | grep -oP 't: \K\d+')
    s=$(echo "$output" | grep -oP 'real\t\K.*' | sed 's/0m//')
    
    echo "$n & $t & $s \\\\"
    echo "--------------------------"
done