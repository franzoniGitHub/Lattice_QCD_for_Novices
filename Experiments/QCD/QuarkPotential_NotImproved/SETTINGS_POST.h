////////////////////////////////////////////////////////////////////////
/// \file SETTINGS_POST.h
/// \brief Header file to set the adjustable parameters for
///        the postprocessing code.
///
/// Set here the parameters of the system, namely: the name of the input file
/// filename, the boolean option to perform a smearing smeared,
/// the number of smearings to apply Nsmearings,
/// the value of the smearing parameter smear_par and
/// the type of analysis my_type.
///
////////////////////////////////////////////////////////////////////////

#ifndef SETTINGS_H
#define SETTINGS_H

#include <string>
#include "source/Metropolis.h"

//*************** PARAMETERS **********************//
std::string filename =
    "DataOutput_8x8x8x8_100NofSU3_10Ncf_notimproved.dat";  ///< Name of the input file
bool smeared = true;                                       ///< Option to perform smearings
int Nsmearings = 1;                                        ///< Number of smearings
double smear_par = 1. / 12.;                               ///< Smearing parameter

/// Type of analysis
///
/// Type of analysis to perform on the input data. The possible choices are the following:
/// Type::PlaquetteRectangule -> Compute 1x1 and 1x2 Wilson loops expectation values
/// Type::QuarkPotential -> Compute quark-quark potential, draw plot and fit using ROOT
/// Type::Custom -> Define in CUSTOM_POST.h a function of the link variables whose expectation value
/// is to be computed \see Type
Type my_type = Type::QuarkPotential;

//************** END PARAMETERS *******************//

#endif
