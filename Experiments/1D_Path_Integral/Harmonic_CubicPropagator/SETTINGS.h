////////////////////////////////////////////////////////////////////////
/// \file SETTINGS.h
/// \brief Header file to set the adjustable parameters of the system
///
/// Set here the parameters of the system, namely: the number of sites
/// N, the number of correlated path configurations to skip
/// Ncorr, the number of path configurations to be sampled
/// Ncf, the number of bootstrap copies NBootStraps,
/// the magnitude of a path update epsilon, the time spacing
/// a and the output filename output_name. If you wish to run the
/// Metropolis::ComputeBinnedEnergyEstimators method, make sure to
/// define also the integer for the bin size.
///
////////////////////////////////////////////////////////////////////////

#ifndef SETTINGS_H
#define SETTINGS_H

//*************** PARAMETERS **********************
int N = 20;             ///< Number of sites
int Ncorr = 20;         ///< Number of correlated path configurations to skip before next acquisition
int Ncf = 10000;        ///< Total number of path configurations to be sampled
int NBootStraps = 100;  ///< Number of bootstraps to perform on the Ncf configurations
double epsilon = 1.4;   ///< Typical magnitude of a path update
double a = 0.5;         ///< Time discretization, i.e. grid spacing
std::string output_name = "Ncf10000_NB100.dat";  ///< Output filename (relative to this directory)
//************** END PARAMETERS *******************

#endif
