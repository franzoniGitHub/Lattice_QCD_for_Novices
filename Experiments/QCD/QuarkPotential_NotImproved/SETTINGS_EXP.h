////////////////////////////////////////////////////////////////////////
/// \file SETTINGS_EXP.h
/// \brief Header file to set the adjustable parameters for
///        the experiment code.
///
/// Set here the parameters of the system, namely: the 4D grid dimension
/// NCells, the grid spacing a (in fm units), the values of beta, beta_tilde
/// and u0 in the Wilson lagrangian and its improved version,
/// the typical magnitude of a path update epsilon,
/// the number of SU3 matrices to be generated NofSU3,
/// the number of correlated path configurations to skip Ncorr,
/// the number of internal updates on each link inner,
/// the number of path configurations to be sampled Ncf,
/// the boolean option improved to select the action to use
/// and the name of the output file filename.
///
////////////////////////////////////////////////////////////////////////

#ifndef SETTINGS_H
#define SETTINGS_H

#include <vector>
#include <string>

//*************** PARAMETERS **********************//
std::vector<int> NCells = {8, 8, 8, 8};  ///< Number of cells in each space-time direction
double a = 0.25;                         ///< Grid spacing in fm units
double beta = 5.5;          ///< Value of the beta=beta_tilde/u_0 in the Wilson lagrangian density
double beta_tilde = 1.719;  ///< Value of the beta_tilde in the Wilson lagrangian density
double u0 = 0.797;          ///< Tadpole improvement: in case of unimproved action (corresponding to
                            ///< improved=false), it doesn't affect calculations
double epsilon = 0.24;      ///< Typical magnitude of a link update
int NofSU3 = 100;           ///< Number of SU3 matrices to be generated and used to update the links
int Ncorr = 50;             ///< Number of correlated configurations to skip before next sampling
int inner = 10;             ///< Number of link updates before moving to the next site
int Ncf = 10;               ///< Number of total configurations to be sampled
bool improved = false;      ///< Do you wish to use the improved wilson action?
std::string filename =
    "DataOutput_8x8x8x8_100NofSU3_10Ncf_notimproved.dat";  ///< Filename of the output data file
//************** END PARAMETERS *******************//

#endif
