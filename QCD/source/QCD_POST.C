////////////////////////////////////////////////////////////////////
///  \file QCD_POST.C
///  \brief Code to analyse the output from QCD_EXP.C and provide
///         Montecarlo estimates
///
///  In the main function, the Metropolis configurations and settings are read
///  from file, the required smearings are applied and statistics are computed.
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
#include "../SETTINGS_POST.h"

using namespace arma;

/// Main function
///
/// \return 0 in a normal excecution, 1 if an exception is caught.
int main()
{
  try {
    auto start = std::chrono::steady_clock::now();

    // Initialize the Metropolis instance by reading from file
    Metropolis latticeQCD(filename);
    // If required, apply the smearing operation
    if (smeared) latticeQCD.SpatialSmearing(Nsmearings, smear_par);
    // Compute the statistics of the required type
    latticeQCD.ComputeStatistics(my_type);

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