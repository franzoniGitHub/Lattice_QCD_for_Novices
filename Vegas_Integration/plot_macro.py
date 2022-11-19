## @file plot_macro.py
# Python code to generate a plot from the
# VegasHarmOsc.C output.
#
# Python code to generate a plot from the
# VegasHarmOsc.C output. This macro also performs
# a fit to the function fitplot_py(), providing an
# estimate for the ground energy E_0. Various graphical
# options are implemented and tuned in the code. Type:\n
# $ python plot_macro.py\n
# to run and generate the plot "plot_py.png".

import re
import numpy as np
import matplotlib.pyplot as plt
from numpy import genfromtxt
from scipy.optimize import least_squares

## Function to read time_bound from header file
#
#   @param filename name of the file to look into.
#   @see SETTINGS.h
def read_param_from_file(filename):
    regexp = re.compile(r' *double time_bound *(=) *([0-9.-]+); *.*')
    time_bound=0.
    with open(filename) as f:
        for line in f:
            match = regexp.match(line)
            if match:
                print("Found in header file time_bound =", match.group(2))
                time_bound=match.group(2)
    return time_bound

## Fit function to be used in the scipy.optimize.least_squares method
#
#   @param par array of parameters, where par[0]=E_0 and par[1]=T
#   @param x coordinate x
#   @param y data to be fitted
#   @see physical_params
def fitfunction_py(par,x,y):
    return (np.exp(-x**2-par[0]*par[1])/np.sqrt(3.14159))-y

## Fit function to be used in the plot with fitted paramaters
#
#   @param x coordinate x
#   @param E_0 ground state energy
#   @param T propagation time
#   @return \f[\left|\left\langle \psi | E_{0}\right\rangle\right|^{2}e^{-E_{0}T}=\frac{1}{\sqrt{\pi}}e^{-x^{2}-E_{0}T}\f]
#   @see physical_params
def fitplot_py(x,E_0,T):
    return (np.exp(-x**2-E_0*T)/np.sqrt(3.14159))

## \cond MAIN
#************** EXECUTION ******************#

input_name="SETTINGS.h"
output_name="plot_py.png"

data=genfromtxt('output_file.dat', delimiter='\t', skip_header=1)
x=data[:,0]
result=data[:,1]
error=data[:,2]
asymptotic=data[:,3]
n_data=x.size

time_bound=read_param_from_file(input_name)

par=np.zeros(2)
par[0]=0.5
par[1]=time_bound
#Fit data to the fitting function using the least square method implemented in the the scipy function
fit_lsq = least_squares(fitfunction_py, par, bounds=([0.,1.],[time_bound,time_bound]), loss='soft_l1', args=(x[:], result[:]))
x_fit = np.linspace(x[0], x[n_data-1],500) #generation of points for the graph
y_fit = fitplot_py(x_fit, *fit_lsq.x) #generation of points for the graph

print("Generating plot with following data")
print("x     result     error      asymptotic")
for i in range(n_data):
    print(x[i],result[i],error[i],asymptotic[i])

# Plot all data and a legend
plt.errorbar(x,result, yerr=error, color='blue', linestyle="None", marker='o', markersize=4, label="Montecarlo data")
plt.plot(x_fit,y_fit, color='green', marker='None', label="Exact asymptotic for $E_{0}=\frac{1}{2}$")
plt.legend(loc='lower center')
# Set title and x-y axis labels for the Axes
plt.title("1D Harmonic Oscillator with Vegas", fontsize=18)
plt.xlabel('$x$')
plt.ylabel('$<x|e^{-HT}|x>$')
# Adds gridlines to the Axes
plt.grid(linestyle='--', color='gray')
# Save the current figure to a file
plt.savefig(output_name)
plt.show()

## \endcond
