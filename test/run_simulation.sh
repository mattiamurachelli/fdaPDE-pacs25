#!/bin/bash

testing=$HOME/fdaPDE-testing
cpp=$testing/fdaPDE-cpp
core=$cpp/fdaPDE/core

cd $testing/test

if [ -f script ]; then rm script; fi
g++ -o script $1 -I$cpp -I$core -I/usr/include/eigen3 -O2 -std=c++20 -march=native -s
./script
