#!/bin/bash

# bash script to compile programs using the gsl library
#****** SETTINGS *******
MAIN="VegasHarmOsc.cpp"
MAIN_DIR="source"
#*** END OF SETTINGS ***

rm executable
cd $MAIN_DIR
rm *.o
g++ -c -Wall -I/usr/local/include -o object.o $MAIN
g++ -L/usr/local/lib object.o -lgsl -lgslcblas -lm -o executable
mv executable ../
cd -
