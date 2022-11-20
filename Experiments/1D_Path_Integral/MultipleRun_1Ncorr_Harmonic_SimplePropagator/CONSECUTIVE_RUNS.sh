#!/bin/bash
g++ source/Metropolis.cpp source/1D_path_integral.cpp -o executable

for (( i = 0 ; i <= 20 ; i += 1 )) ; do
  echo "Step n.$i"
  ./executable > waste.log
  python plot_macro.py > waste.log
  mv plot_py.png plot_py_$i.png
  mv output.dat output_$i.dat
done
