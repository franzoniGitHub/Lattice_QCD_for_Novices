////////////////////////////////////////////////////////////////////
///  \file VegasHarmOsc.C
///  \brief Code to path integrate a harmonic oscillator by direct integration
///
///  In the main function, the physical parameters are set and the
///  loop path integration is performed using the GSL implementation
///  of the Vega routine. The integration is repeated using different
///  values of the start-end point x, so as to provide data to approximate
///  the shape of the ground state function of the system.
///  See also the GSL manual at https://www.gnu.org/software/gsl/doc/html/index.html.
///  Follow the comments in the source code of the main function below for a more detailed
///  description. \include MainForDoxygen.C
///
////////////////////////////////////////////////////////////////////

#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <fstream>
#include <iomanip>
#include <iostream>

////////////////////////////////////////////////////////////////////
///  Parameters of the physical system
///
///  Structure containing the parameters of the system.
///
////////////////////////////////////////////////////////////////////
struct physical_params
{
  double x;     ///< starting and ending position of the loop path
  double T;     ///< propagation time
  double mass;  ///< mass of the particle
};

////////////////////////////////////////////////////////////////////
///  Potential
///
///  Function defining the shape of the potential felt by the
///  particle: modify this to change the interaction.
///  \param x position
///  \param fp pointer to the set of parameters
///  \return value of the potential computed in \a x
///  \see physical_params
///
////////////////////////////////////////////////////////////////////
double V(double x, physical_params* fp)
{
  return pow(x, 2.0) * (fp->mass) / 2.0;
}

////////////////////////////////////////////////////////////////////
///  Integrand of the path integral
///
///  Function to compute the integrand in the path
///  integral. When evaluating the action, periodic boundary conditions
///  are implemented. The signature of this function is dictated by its
///  use within the GSL library.
///  \param path array containing the discretized 1D path
///  \param dim {size of the path array (due to periodic boundary
///             conditions, it corresponds to N-1, where N is the
///             lattice size)}
///  \param params pointer to the set of physical parameters
///  \return \f$Ae^{-S}\f$, where \f$S\f$ is the action evaluated on the path
///          and \f$A=(\frac{m}{2\pi a})^{N/2}\f$, \f$m\f$ being the mass and
///          \f$a=\frac{T}{N}\f$.
///  \see physical_params
///
////////////////////////////////////////////////////////////////////
double integrand(double* path, size_t dim, void* params)
{
  struct physical_params* gp = (struct physical_params*)params;

  int N = (int)dim + 1;
  double a = (gp->T) / ((double)N);
  double m = gp->mass;
  double A = pow(m / (2.0 * M_PI * a), (double)N / 2.0);

  double action = 0.0;
  for (int i = 0; i < N; i++) {
    if (i == 0) action += (m * pow(path[0] - gp->x, 2.0) / (2 * a)) + a * V(gp->x, gp);
    if (i > 0 && i < N - 1)
      action += (m * pow(path[i] - path[i - 1], 2.0) / (2 * a)) + a * V(path[i - 1], gp);
    if (i == N - 1)
      action += (m * pow(gp->x - path[N - 2], 2.0) / (2 * a)) + a * V(path[N - 2], gp);
  }
  return A * exp(-action);
}

////////////////////////////////////////////////////////////////////
///  Asymptotic exact function
///
///  Asymptotic exact function in the large time limit computed from standard quantum mechanics
///  \param x position
///  \param fp pointer to the set of parameters
///  \return \f[\left|\left\langle \psi |
///  E_{0}\right\rangle\right|^{2}e^{-E_{0}T}=\frac{1}{\sqrt{\pi}}e^{-x^{2}-E_{0}T}\f] \see
///  physical_params
///
////////////////////////////////////////////////////////////////////
double asymptotic(double x, physical_params* fp)
{
  return exp(-x * x - 0.5 * (fp->T)) / sqrt(3.14159);
}

int main()
{
#include "../SETTINGS.h"  // Include the parameters definition

  double res, err;       // result and error outputs
  double xl[N_dim - 1];  // lower bound vector
  double xu[N_dim - 1];  // upper bound vector

  // initialization of the bound vectors
  for (int i = 0; i < N_dim - 1; i++) {
    xl[i] = -space_bound;
    xu[i] = space_bound;
  }
  size_t dimension = (size_t)(N_dim - 1);

  // Definition and initialization of GSL random generation variables
  // These are necessary for the Vegas routine, which uses Montecarlo importance sampling
  const gsl_rng_type* T;
  gsl_rng* my_gsl_rng;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  my_gsl_rng = gsl_rng_alloc(T);

  // Prepare for the output on standard output with a header
  std::cout << std::fixed << std::setprecision(8) << "x\tpath integral\tsigma\texact asymptotic\n";

  // Open an output file and print a header
  std::ofstream out_file("output_file.dat");
  out_file << std::fixed << std::setprecision(8) << "x\tpath integral\tsigma\texact asymptotic\n";

  // Evaluate the path integral for different cyclic boundary conditions contained in my_fp
  double total_integral = 0.0;
  for (double x_now = -x_loop_bound; x_now <= x_loop_bound + 1.e-8; x_now += x_loop_step) {
    size_t calls = 500000;  // Number of calls of the integration routine
    struct physical_params my_fp = {x_now, time_bound, mass};  // Setting the physical parameters

    // Definition of the function with parameters for Montecarlo integration
    gsl_monte_function my_gsl_integrand = {&integrand, dimension, &my_fp};

    // Definition and initialization of a workspace to maintain the state of integration
    gsl_monte_vegas_state* my_gsl_state = gsl_monte_vegas_alloc(dimension);

    // Vegas routine to integrate over the hypercubic region whose lower and upper limits are
    // defined by xl and xup, respectively Do 10000 iterations and store result and error in res and
    // error
    gsl_monte_vegas_integrate(
        &my_gsl_integrand, xl, xu, dimension, 10000, my_gsl_rng, my_gsl_state, &res, &err);
    do  // Iterate the integration until the required accuracy target is met, namely, when
        // Chisquare/#dofs is sufficiently close to 1
    {
      gsl_monte_vegas_integrate(
          &my_gsl_integrand, xl, xu, dimension, calls / 5, my_gsl_rng, my_gsl_state, &res, &err);
    } while (fabs(gsl_monte_vegas_chisq(my_gsl_state) - 1.0) > 0.5);

    // Print all the results for x=x_now
    std::cout << x_now << "\t" << res << "\t" << err << "\t" << asymptotic(x_now, &my_fp)
              << std::endl;
    out_file << x_now << "\t" << res << "\t" << err << "\t" << asymptotic(x_now, &my_fp)
             << std::endl;

    // Update a rough estimate of the integral of the path integral over the start-end points x
    total_integral += x_loop_step * res;

    // Free the memory allocated
    gsl_monte_vegas_free(my_gsl_state);
  }
  std::cout << "Estimate of the zero energy level E_0=" << -(log(total_integral)) / time_bound
            << std::endl;

  // Free the memory allocated and close the output file
  gsl_rng_free(my_gsl_rng);
  out_file.close();
  return 0;
}

////////////////////////////////////////////////////////////////////
/// \mainpage
///
/// \section intro Introduction
/// The purpose of this code is to reproduce the first excercise in Lepage's
/// article "Lattice QCD for Novices". A harmonic oscillator is explicitely
/// path integrated in order to provide estimates of its zero energy and eigenfunction.
/// The integration is performed via the Vega routine, as implemented in the GSL
/// library.
///
/// The source code is composed of two files: SETTINGS.h and VegasHarmOsc.C. \n
/// The plot routines are given in two versions: plot_macro.C and plot_macro.py.
///
/// \section reqs Requirements
/// - To compute the results: g++ compiler and working GSL installation (download at
/// https://www.gnu.org/software/gsl/)
/// - To plot the results: either ROOT (download at https://root.cern/install/) or a working Python
/// installation.
///
/// \section build How to Build and Run
///  - Check and set the parameters in SETTINGS.h
///  - Execute the script BUILD.sh\n
///    $ ./BUILD.sh
///  - Run the executable\n
///    $ ./executable
///  - Check the output file
///
/// \section plot How to plot
///  - Using ROOT\n
///    $ root plot_macro.C
///  - Using Python\n
///    $ python plot_macro.py
///
/// \section clean How to Clean
///  - Remove all the output \n
///    $ ./CLEAN.sh
///  - NOTE: plots and output files are deleted in this procedure
///
////////////////////////////////////////////////////////////////////
