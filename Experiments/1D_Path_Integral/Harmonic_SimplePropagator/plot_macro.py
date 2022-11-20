## @file plot_macro.py
#Python code to generate a plot from the
# 1D_path_integral.cpp output.
#
#Python code to generate a plot from the
# 1D_path_integral.cpp output of the method
#Metropolis::ComputeEnergyEstimators().Various graphical
#options are implemented and tuned in the code.Type :\n
#$ python plot_macro.py\n
#to run and generate the plot "plot_py.png".

import re
import numpy as np
import matplotlib.pyplot as plt
from numpy import genfromtxt
from scipy.optimize import least_squares

## Function to read parameters a, output_name and N from header file
#
#@param filename name of the header file to look into.
#@see SETTINGS.h
def read_param_from_file(filename):
    regexp_N = re.compile(r' *int N *(=) *([0-9.-]+) *; *.*')
    regexp_a = re.compile(r' *double a *(=) *([0-9.-]+) *; *.*')
    regexp_output_name = re.compile(r' *std::string output_name *(=) *"(\w*[.]\w*)" *; *.*')
    N=0
    a=0.
    output_name=""
    f=open(filename)
    for line in f:
        match_N = regexp_N.match(line)
        if match_N:
            print("Found in header file N =", match_N.group(2))
            N=match_N.group(2)
    f=open(filename)
    for line in f:
        match_a = regexp_a.match(line)
        if match_a:
            print("Found in header file a =", match_a.group(2))
            a=match_a.group(2)
    f=open(filename)
    for line in f:
        match_output_name = regexp_output_name.match(line)
        if match_output_name:
            print("Found in header file output_name =", match_output_name.group(2))
            output_name=match_output_name.group(2)
    return N,a,output_name

## \cond MAIN
#* * * * * * * * * * * * * * EXECUTION * * * * * * * * * * * * * * * * * * #

input_name="SETTINGS.h"
plot_name="plot_py.png"

N,a,output_name=read_param_from_file(input_name)
a=float(a)

data=genfromtxt(output_name, delimiter='\t', skip_header=1)
n_data=round(0.34*(data[:,0].size))+1
time=data[:n_data,0]
result=data[:n_data,3]
error=data[:n_data,4]

print("Generating plot with following data")
print("time     deltaE     error")
for i in range(n_data):
    print(time[i],result[i],error[i])

#Generate the exact constant line deltaE = 1
dE_exact=[1.,1.]
t_exact=[-0.5,n_data*a*1.5]

#Plot all data and a legend
plt.errorbar(time,result, yerr=error, color='blue', linestyle="None", marker='o', markersize=4, label="Montecarlo data")
plt.plot(t_exact,dE_exact, color='red', marker='None', label="Exact asymptotic $\Delta E=1$")
plt.legend(loc='upper right')
#Set title and x - y axis labels for the Axes
plt.title("1D Harmonic Oscillator, Simple Propagator", fontsize=18)
plt.axis([-0.2,(n_data-1)*a+0.2,0.,2.])
plt.xlabel('$t$')
plt.ylabel('$\Delta E(t)$')
#Adds gridlines to the Axes
plt.grid(linestyle='--', color='gray')
#Save the current figure to a file
plt.savefig(plot_name)
plt.show()

## \endcond
