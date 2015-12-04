#!/bin/sh
make
./main 1> output 2> residual
./plot.py
