////////////////////////////////////////////////////////////////////
///  \file QCD_EXP.C
///  \brief Code to path integrate a gluonic system with Montecarlo method
///
///  In the main function, the physical parameters are set and the
///  path integration is performed using the Metropolis algorithm,
///  as implemented in the Metropolis class.
///  Then, the ensemble of configurations is printed on file.
///  Follow the comments in the source code of the Metropolis class for
///  a more detailed description.
///
////////////////////////////////////////////////////////////////////
#include <armadillo>
#include <chrono>
#include <iostream>
#include <string>
#include <vector>
#include "my4Vector.h"
#include "Path.h"
#include "Metropolis.h"
#include "../SETTINGS_EXP.h"
using namespace arma;

/// Main function
///
/// \return 0 in a normal excecution, 1 if an exception is caught.
int main()
{
  try {
    auto start = std::chrono::steady_clock::now();

    // Define the parameter container for the Metropolis constructor
    std::vector<int> int_params = {NofSU3, Ncorr, inner, Ncf};
    std::vector<double> double_params = {a, beta, beta_tilde, u0, epsilon};

    // Initialize the Metropolis instance
    Metropolis latticeQCD(NCells, int_params, double_params, improved);
    // Run the Metropolis algorithm to generate physical configurations
    latticeQCD.RunMetropolis();
    // Print the settings and the path configurations on file
    latticeQCD.PrintAllOnFile(filename);

    // Print the execution time
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    double deltaT_sec = elapsed_seconds.count();
    std::cout << "Execution time: " << deltaT_sec / 60. << " min\n";

    return 0;
  } catch (...) {
    std::cout << "Exception caught.\n";
    return 1;
  }
}

////////////////////////////////////////////////////////////////////
/// \mainpage
///
/// \section intro Introduction
/// The purpose of this code is to reproduce the third set of excercises in Lepage's
/// article "Lattice QCD for Novices". In the simulation phase, gluonic path integrals are evaluated
/// with the Metropolis algorithm, using two possible versions of the Wilson action: the standard
/// Wilson action and its improved version. The standard Wilson action is defined as follows:
/// \f[S=-\beta\sum_{x,\mu>\nu}{P_{\mu\nu}}\f]
/// where \f$P_{\mu\nu}\f$ is the square 1x1 Wilson loop,
/// \f$\beta=\tilde{\beta}/u_{0}^{4}\f$ and \f$u_{0}^{4}\f$ is the mean link value.
/// Denoting with \f$R_{\mu\nu}\f$ the 1x2 rectangular Wilson loops, the improved
/// Wilson action reads as follows:
/// \f[S=-\tilde{\beta}\sum_{x,\mu>\nu}{\frac{5}{3}\frac{P_{\mu\nu}}{u_{0}^{4}}-
///                                     \frac{R_{\mu\nu}+R_{\nu\mu}}{12u_{0}^{6}}}\f]
/// For the analysis phase, three options are available for statistics
/// computations: the evaluation of 1x1 and 1x2 Wilson loops, the computation of the quark-quark
/// potential or the expectation value of a user-defined function of the links (as defined
/// in CUSTOM_POST.h). The possibility to split the simulation and the analysis
/// phases is useful to reuse the same Metropolis ensemble for subsequent computations
/// of different observables, thus making the process less time consuming.
///
/// The source code is composed of two programs, QCD_EXP.C and QCD_POST.C. Their purpose is,
/// respectively, to perform an experiment and to later analyse the results. Both these codes are
/// based on the classes my4Vector (for the implementation of operations on position 4-vectors which
/// include periodic boundary conditions), Path (for the definition of a generic lattice
/// configuration) and Metropolis (for the definition of a Metropolis algorithm working on a 4D
/// quantum system with SU3 gauge-symmetry). The classes Path and Metropolis, both rely on the
/// Armadillo library for the use of complex matrices and operations on them. This library is
/// available under the Apache licence (https://opensource.org/licenses/Apache-2.0) and can be
/// downloaded at https://arma.sourceforge.net/download.html. The Armadillo library
/// (http://dx.doi.org/10.21105/joss.00026, https://doi.org/10.1007/978-3-319-96418-8_50) was
/// adopted for its performances, in particular for its optimisation of the matrix operations and
/// the reduction of temporaries.
///
/// The plot routines are given in two versions: plot_macro.C and plot_macro.py.
///
/// \section reqs Requirements
/// - To compute the results: g++ compiler and a working installation of the Armadillo library
/// (download at https://arma.sourceforge.net/download.html)
/// - To plot the results: either ROOT (download at https://root.cern/install/) or a working Python
/// installation.
///
/// \section build How to Build and Run
///  Experiment phase:
///  - It is suggested to create a copy of the QCD directory, so as to keep
///    the sequence of experiments separated
///  - Check and set the parameters in SETTINGS_EXP.h
///  - Execute the script BUILD_EXP.sh \n
///    $ ./BUILD_EXP.sh
///  - Run the executable \n
///    $ ./executable_EXP
///  - Check the output file
///
///  Postprocessing phase:
///  - Check and set the parameters in SETTINGS_POST.h
///  - Execute the script BUILD_POST.sh \n
///    $ ./BUILD_POST.sh
///  - Run the executable \n
///    $ ./executable_POST
///  - Check the output files
///
/// \section plot How to plot
///  - Using ROOT \n
///    $ root plot_macro.C
///  - Using Python \n
///    $ python plot_macro.py
///
/// \section clean How to Clean
///  - Remove all the executables \n
///    $ ./CLEAN.sh
///  - NOTE: plots and output files are NOT deleted in this procedure
///
////////////////////////////////////////////////////////////////////