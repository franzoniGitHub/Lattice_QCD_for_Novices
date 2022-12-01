## @file plot_macro.py
#
#  @brief Python macro to generate a plot from the
#         QCD_POST.cpp output.
#
#  Python macro to generate a plot from the
#  QCD_POST.cpp output of the method
#  Metropolis::ComputeRxTWilsonLoops which is stored in the file
#  RXT_potential_plot_file.dat. Various graphical
#  options are implemented and tuned in the code. Set the maximum
#  value of r/a (r_max) to be displayed in the PARAMETER section
#  of the script. Type:\n
#  $ python plot_macro.py\n
#  to run and generate the plot "QuarkPotential_py.png".
import numpy as np
import matplotlib.pyplot as plt
from numpy import genfromtxt
from scipy.optimize import least_squares

## Fit function to be used in the scipy.optimize.least_squares method
#
#   @param par array of parameters, where par[0]=\f$\sigma\f$, par[1]=b, par[2]=c.
#   @param x coordinate x
#   @param y data to be fitted
def fitfunction_py(par,x,y):
    return par[0]*x-par[1]/x+par[2]-y

## Fit function to be used in the plot with fitted paramaters
#
#   @param r coordinate r
#   @param sigma parameter \f$\sigma\f$
#   @param b parameter b
#   @param c parameter c
#   @return \f[V(r)=\sigma r-\frac{b}{r}+c\f]
def fitplot_py(r,sigma,b,c):
    return sigma*r-b/r+c

## \cond MAIN
#************** EXECUTION ******************#

#************** PARAMETER ******************#
r_max=3
#*******************************************#

input_name="SETTINGS.h"
output_name="QuarkPotential_py.png"

data=genfromtxt('RXT_potential_plot_file.dat', delimiter='\t', skip_header=15)
x=data[:r_max,0]
result=data[:r_max,1]
error=data[:r_max,2]
n_data=x.size

par=np.zeros(3)
par[0]=5.
par[1]=5.
par[2]=0.
#Fit data to the fitting function using the least square method implemented in the the scipy function
fit_lsq = least_squares(fitfunction_py, par, loss='soft_l1', args=(x[:], result[:])) #bounds=([-10.,10.],[-10.,10.],[-10.,10.]),
x_fit = np.linspace(x[0], x[n_data-1],500) #generation of points for the graph
y_fit = fitplot_py(x_fit, *fit_lsq.x) #generation of points for the graph

print("Generating plot with following data")
print("x     result     error")
for i in range(n_data):
    print(x[i],result[i],error[i])

# Plot all data and a legend
plt.errorbar(x,result, yerr=error, color='blue', linestyle="None", marker='o', markersize=4, label="Montecarlo data")
plt.plot(x_fit,y_fit, color='green', marker='None', label="Fit to $V(r)=\sigma r-b/r+c$")
plt.legend(loc='upper center')
# Set title and x-y axis labels for the Axes
plt.title("Quark Potential", fontsize=18)
plt.xlabel('$r/a$')
plt.ylabel('$aV(r)$')
# Adds gridlines to the Axes
plt.grid(linestyle='--', color='gray')
# Save the current figure to a file
plt.savefig(output_name)
plt.show()

## \endcond