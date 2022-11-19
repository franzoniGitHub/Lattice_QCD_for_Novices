////////////////////////////////////////////////////////////////////////
/// \file SETTINGS.h
/// \brief Header file to set the adjustable parameters of the system
///
/// Set here the parameters of the system, namely: the number of sites
/// (\link N \endlink), the number of correlated path configurations to skip
/// (\link Ncorr \endlink), the number of path configurations to be sampled
/// (\link Ncf \endlink), the number of bootstrap copies (\link NBootStraps \endlink)
/// the magnitude of a path update (\link epsilon \endlink), the time spacing
/// (\link a \endlink) and the output filename (\link output_name \endlink)
///
////////////////////////////////////////////////////////////////////////

#ifndef SETTINGS_H
#define SETTINGS_H

//*************** PARAMETERS **********************
int N = 20;            ///< Number of sites
int Ncorr = 1;         ///< Number of correlated path configurations to skip before next acquisition
int Ncf = 1000;        ///< Total number of path configurations to be sampled
int NBootStraps = 50;  ///< Number of bootstraps to perform on the Ncf configurations
int bin_width = 20;    ///< Bin width
double epsilon = 1.4;  ///< Typical magnitude of a path update
double a = 0.5;        ///< Time discretization, i.e. grid spacing
std::string output_name = "output.dat";  ///< Output filename (relative to this directory)
//************** END PARAMETERS *******************

#endif
