////////////////////////////////////////////////////////////////////
///  \file 1D_path_integral.C
///  \brief Code to path integrate a 1D quantum system with Montecarlo method
///
///  In the main function, the physical parameters are set and the
///  loop path integration is performed using the Metropolis algorithm,
///  as implemented in the Metropolis class.
///  Then, the requested observables are computed on the ensemble.
///  Follow the comments in the source code of the Metropolis class for
///  a more detailed description.
///
////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>
#include "Metropolis.h"

int main()
{
  try {
#include "../SETTINGS.h"

    // Initialize the Metropolis instance
    Metropolis myAlgorithm(N, Ncorr, Ncf, NBootStraps, epsilon, a);
    // Run the Metropolis algorithm to generate physical configurations
    myAlgorithm.RunMetropolis();
    // Compute both the gamma and energy estimators, print and save them
    myAlgorithm.ComputeEnergyEstimators(output_name);

    return 0;
  } catch (...) {
    std::cout << "Exception caught\n";
    return 1;
  }
}

////////////////////////////////////////////////////////////////////
/// \mainpage
///
/// \section intro Introduction
/// The purpose of this code is to reproduce the second excercise in Lepage's
/// article "Lattice QCD for Novices". A 1D quantum system is path integrated
/// using the Montecarlo algorithm in order to provide estimates of first
/// excited state energy gap.
///
/// The source code is composed of five files: SETTINGS.h,
/// CUSTOM.h, Metropolis.h, Metropolis.C and
/// 1D_path_integral.C \n
/// The plot routines are given in two versions: plot_macro.C and plot_macro.py.
///
/// \section reqs Requirements
/// - To compute the results: g++ compiler
/// - To plot the results: either ROOT (download at https://root.cern/install/) or a working Python
/// installation.
///
/// \section build How to Build and Run
///  - It is suggested to create a copy of the 1D_path_integral directory, so as to keep
///    the sequence of experiments separated
///  - Provide the system definition and the required observables in CUSTOM.h
///  - Check and set the parameters in SETTINGS.h
///  - Execute the script BUILD.sh \n
///    $ ./BUILD.sh
///  - Run the executable \n
///    $ ./executable
///  - Check the output file
///
/// \section plot How to plot
///  - Using ROOT \n
///    $ root plot_macro.C
///  - Using Python \n
///    $ python plot_macro.py
///
/// \section clean How to Clean
///  - Remove all the output \n
///    $ ./CLEAN.sh
///  - NOTE: plots and output files are deleted in this procedure
///
////////////////////////////////////////////////////////////////////
