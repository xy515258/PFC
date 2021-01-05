#!/bin/sh
for step in {0..100}
do
    let a=$step*5000
    echo $a
    g++ -std=c++11 find_peak.cpp
    ./a.out "$a"
done
