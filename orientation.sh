#!/bin/sh
module unload rh/devtoolset/9 

for step in {0..100}
do
    let a=340000+$step*20000
    echo $a
    g++ -std=c++11 orientation_analysis.cpp
    ./a.out "$a"
done
