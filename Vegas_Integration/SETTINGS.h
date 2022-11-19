////////////////////////////////////////////////////////////////////////
/// \file SETTINGS.h
/// \brief Header file to set the adjustable parameters of the system
///
/// Set here the parameters of the system, namely: the space boundaries
/// (\link space_bound \endlink), the propagation time (\link time_bound \endlink), the
/// lattice dimension (\link N_dim \endlink) and the particle mass (\link mass \endlink).
///
////////////////////////////////////////////////////////////////////////

#ifndef PARAMS_H
#define PARAMS_H

//***** PARAMETERS **********************************************
double space_bound = 5;     ///< the path integration is limited on the interval [-space_bound,
                            ///< space_bound] to save computational time
double time_bound = 4;      ///< propagation time of the particle
int N_dim = 8;              ///< dimension of the 1D discretized path
double x_loop_bound = 3.0;  ///< The start-end points of the loop paths are sampled in the interval
                            ///< [-x_loop_bound, x_loop_bound]
double x_loop_step = 0.2;   ///< Step between two successive start-end points of the loop paths
                            ///< sampled in [-x_loop_bound, x_loop_bound]
double mass = 1.0;          ///< particle mass
//***** END OF PARAMETERS ***************************************

#endif
