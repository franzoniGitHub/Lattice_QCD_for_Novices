////////////////////////////////////////////////////////////////////////
/// \file Metropolis.h
/// \brief Header file for the definition of the class Metropolis
///
/// Header file containing the definitions of the attributes and members
/// of the class Metropolis. Further comments may be found in the
/// implementation file of this class.
////////////////////////////////////////////////////////////////////////
#ifndef METROPOLIS_H
#define METROPOLIS_H

#include <armadillo>
#include <complex>
#include <random>
#include <string>
#include <vector>
#include "Path.h"
using namespace arma;

/// Type enum class
///
/// Enum class which collects the available types of analysis for the method
/// Metropolis::ComputeStatistics.
enum class Type {
  PlaquetteRectangule,  ///< Select the 1x1 and 1x2 Wilson loops expectation values
  QuarkPotential,       ///< Select the quark potential computation
  Custom                ///< Select a customized analysis, as defined in CUSTOM_POST.h
};

/// Metropolis class
///
/// The instances of this class represent an
/// implementation of the Metropolis algorithm for a
/// SU(3) gauge-symmetric physical quantum system in 4D. The physical system is
/// described by its euclidean action S, defined in the function
/// Metropolis::S in terms of link variables. In the current implementation, the
/// available actions are the standard Wilson action and its improved version.
/// The Metropolis algorithm is run with the Metropolis::RunMetropolis
/// method. The configurations obtained can be analysed by choosing a Type in the
/// Metropolis::ComputeStatistics method, or printed on file with
/// Metropolis::PrintAllOnFile for later postprocessing and analysis.
/// \see intro
class Metropolis
{
 private:
  int fNofSU3;       ///< Number of SU3 matrices to be generated and used to update the links
  int fNcorr;        ///< Number of correlated configurations to skip before next sampling
  int fInnerCycles;  ///< Number of link updates before moving to the next site
  int fNcf;          ///< Number of total configurations to be sampled

  double fA;          ///< Grid spacing in fm units
  double fBeta;       ///< Value of the beta=beta_tilde/u_0 in the Wilson lagrangian density
  double fBetaTilde;  ///< Value of the beta_tilde in the Wilson lagrangian density
  double fU0;         ///< Tadpole improvement: in case of unimproved action (corresponding to
                      ///< improved=false), fBeta is used, and this variable doesn't affect
					  ///< calculations
  double fEpsilon;    ///< Typical magnitude of a link update

  bool fImproved;  ///< Boolean option to select the action to use, either the improved (true) or
                   ///< standard (false) Wilson action

  std::vector<cx_dmat> fSetOfSU3;  ///< Set of SU3 matrices used to update the links
  Path fPath;                      ///< Path object which defines the current lattice configuration
  std::vector<Path> fResult;  ///< Vector of Path configurations: it stores the Montecarlo ensemble
                              ///< obtained by running Metropolis::RunMetropolis

  /************************ Private Methods ***************************/
  /// Gamma
  ///
  /// Auxiliary function to compute the terms of the standard Wilson action variation independent
  /// on the updated link.
  /// \param x 4-position where Gamma is evaluated
  /// \param mu value of the polarization on which Gamma is evaluated
  cx_dmat Gamma(const my4Vector& x, int mu) const;

  /// Gamma improved
  ///
  /// Auxiliary function to compute the terms of the improved Wilson action variation independent
  /// on the updated link.
  /// \param x 4-position where the improved Gamma is evaluated
  /// \param mu value of the polarization on which the improved Gamma is evaluated
  cx_dmat GammaImproved(const my4Vector& x, int mu) const;

  /// Action
  ///
  /// Action S at point x due to \f$U_{\mu}(x)\f$; the method distinguishes the improved
  /// and the unimproved action depending on the value of fImproved.
  /// \param x 4-position where the action terms are evaluated
  /// \param mu value of the polarization on which the action terms evaluated
  /// \param gamma_x_mu output of the method Metropolis::Gamma at x and mu.
  /// \param gamma_improved_x_mu output of the method Metropolis::GammaImproved at x and mu. If
  /// fImproved=false, this entry is unused. \see Metropolis::Gamma, Metropolis::GammaImproved, \ref intro
  double S(const my4Vector& x,
           int mu,
           const cx_dmat& gamma_x_mu,
           const cx_dmat& gamma_improved_x_mu) const;

  /// Print status
  ///
  /// Auxiliary method to print a status bar when performing time demanding loops
  /// \param index loop index of the current iteration
  /// \param max_index max index value of the loop
  void PrintStatus(int index, int max_index) const;

  /// Print settings on file
  ///
  /// Auxiliary method to print the experiment settings as a header of an output file.
  /// \param file output file object on which the settings are printed
  void PrintSettingsOnFile(std::ofstream& file) const;

 public:
  Metropolis() = delete;

  /// Standard Metropolis constructor from SETTINGS_EXP.h parameters
  ///
  /// \param N vector for the grid dimensions in the 4D.
  /// \param integer_params vector for the integer parameters in the form {fNofSU3, fNcorr,
  /// fInnerCycles, fNcf} \param floating_params vector for the floating parameters in the form {fA,
  /// fBeta, fBetaTilde, fU0, fEpsilon} \param isimproved boolean option to select the action to
  /// use, either the improved (true) or standard (false) Wilson action
  Metropolis(std::vector<int> N,
             std::vector<int> integer_params,
             std::vector<double> floating_params,
             bool isimproved);

  /// Metropolis constructor from file
  ///
  /// \param infile input file: it should be a previous output of the method
  /// Metropolis::PrintAllOnFile in the non-verbose option \see Metropolis::PrintAllOnFile
  Metropolis(std::string infile);

  /// Metropolis destructor
  ~Metropolis();

  /// \return vector containing the lattice dimensions in the 4D
  std::vector<int> GetNCells() const;

  /// \return fNcorr
  int GetNcorr() const;

  /// \return fNcf
  int GetNcf() const;

  /// \return fEpsilon
  double GetEpsilon() const;

  /// \return fBeta
  double GetBeta() const;

  /// \return fBetaTilde
  double GetBetaTilde() const;

  /// \return fU0
  double GetU0() const;

  /// \return fImproved
  bool IsImproved() const;

  /// \return fSetOfSU3
  std::vector<cx_dmat> GetSetOfSU3() const;

  /// \return fPath, the current lattice configuration
  Path GetCurrentPath() const;

  /// \return fResult, the vector containing the Metropolis ensemble of lattice configurations
  std::vector<Path> GetCurrentResult() const;

  /// Print current path on standard output
  void PrintPathOnScreen() const;

  /// Print current path on file
  /// \param outfile name of the file to print on
  void PrintPathOnFile(std::string outfile) const;

  /// Print all on file
  ///
  /// Print settings and Montecarlo configurations on file.
  /// \param outfile name of the file to print on
  /// \param opt option for the type of output:\n
  ///           - by default, the print mode is silent, which is the required option for later
  ///           postprocessing (matrix elements are printed with maximum accuracy, no headers)\n
  ///           - as an option, 'V' may be specified for a verbose print, which is more user
  ///           readable and it saves disk space, though it is less accurate
  void PrintAllOnFile(std::string outfile, char opt = 'S') const;

  /// Clear the result
  ///
  /// Bring the result vector fResult to its default value
  void ClearResult();

  /// Update the current fPath configuration
  /// \return acceptance ratio for the update on the lattice
  double UpdateCurrentPath();

  /// Spatial smearing
  ///
  /// Perform Ntimes a spatial smearing using a discretized version of the
  /// gauge-covariant derivative on the set of Montecarlo configurations in fResult.
  /// \param Ntimes number of consecutive spatial smearings to apply on the lattice
  /// \param smearing_par value of the smearing parameter
  /// \see Metropolis::GaugeDerivative
  void SpatialSmearing(int Ntimes, double smearing_par);

  /// Compute the discretized gauge derivative
  ///
  /// Compute the discretized gauge-covariant derivative summed over all directions on the i-th
  /// fResult element. \param i index of the Metropolis configuration in fResult the derivation is
  /// applied on \return object of type Path containing the gauge derivative on the lattice
  /// configuration
  Path GaugeDerivative(int i) const;

  /// Run the Metropolis algorithm
  ///
  /// Perform a complete run of the Metropolis algorithm and collect the set of configurations in
  /// the fResult vector
  void RunMetropolis();

  /// Randomize the SU3 matrices
  ///
  /// Fill the fSetOfSU3 set with a set of SU3 random matrices
  void RandomizeSU3();

  /// Evaluate a Wilson loop on a Metropolis configuration
  ///
  /// Compute a N_mu X N_nu Wilson loop in the plane mu - nu  at position x using the j-th
  /// configuration in fResult. \param N_mu length of the loop in the mu direction in lattice units
  /// (integer) \param N_nu length of the loop in the nu direction in lattice units (integer) \param
  /// mu direction of one of the sides of the loop \param nu direction of one of the sides of the
  /// loop \param x position of one of the corners of the loop \param j index of the j-th path
  /// configuration on which the loop is evaluated
  double WilsonLoop(int N_mu, int N_nu, int mu, int nu, const my4Vector& x, int j) const;

  /// Compute the 1x1 and 1x2 Wilson loop expectation values
  ///
  /// Compute plaquette (1x1) and rectangule (1x2) expectation values on the current fResult
  /// and print the results on the standard output.
  void ComputePlaquetteRectangule() const;

  /// Compute the RxT Wilson loop expectation values
  ///
  /// Compute the RxT Wilson loop expectation values in planes where one is a spatial direction and
  /// one is the time direction. The results are then used to evaluate the ratio W(R,T)/W(R,T+fA)
  /// and an estimate of the quark-quark potential using the time-asymptotic value of the ratio. All
  /// the quantities are printed in separate files, respectively, RXT_loops_file.dat,
  /// RXT_potential_file.dat and RXT_potential_plot_file.dat. The last file is used as input file in
  /// the ROOT and python plotting macros. \see plot_macro.C, plot_macro.py
  void ComputeRxTWilsonLoops() const;

  /// Compute the expectation value of a custom function
  ///
  /// Compute the expectation value of a custom function. The definition of the function must be
  /// given in the CUSTOM_POST.h file following the indications in the comments.
  void ComputeCustom() const;

  /// Compute statistics of a predefined type
  ///
  /// Compute the statistical analysis of a predefined type in the enum class Type
  /// \param mytype type of statistical analysis, as defined in the enum class Type
  void ComputeStatistics(Type mytype) const;
};

#endif
