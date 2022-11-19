////////////////////////////////////////////////////////////////////////
/// \file Metropolis.h
/// \brief Header file for the definition of the class Metropolis
///
/// Header file containing the definitions of the attributes and members
/// of the class Metropolis. The instances of this class represent an
/// implementation of the Metropolis algorithm for a physical quantum
/// system in 1D via path integral formalism. The physical system is
/// described by its euclidean action S, defined in the function
/// Metropolis::S and implemented in the CUSTOM.h file.
/// The observable Gamma, whose implementation is in the same CUSTOM.h
/// file, is estimated by averaging over the Montecarlo ensemble:
/// the estimators for the chosen Gamma are readily provided by the method
/// Metropolis::ComputeEstimators.
/// In addition, an estimate of the energy gap as a function of time can be
/// obtained via Metropolis::ComputeEnergyEstimators using the propagator
/// defined in CUSTOM.h: in this case, averages and standard
/// deviations of the bootstrap copies of the Montecarlo ensemble are used.
/// These separate implementations aim at providing flexibility both for
/// choice of the physical system and that of the observable.
/// Further comments may be found in the implementation file of this class.
////////////////////////////////////////////////////////////////////////

#ifndef METROPOLIS_H
#define METROPOLIS_H

#include <vector>

//////////////////////////////////////////////////////////////////////////
/// Metropolis class
///
/// Class representing a 1D quantum physical system whose evolution can be
/// computed by path integration using the Metropolis algorithm.
//////////////////////////////////////////////////////////////////////////
class Metropolis
{
 private:
  const int fN;            ///< Number of grid points for the discretized trajectory
  const int fNcorr;        ///< Number of correlated configurations to skip before next acquisition
  const int fNcf;          ///< Number of configurations to sample
  const int fNBootStraps;  ///< Number of statistical bootstraps to perform
  const double fEpsilon;   ///< Typical magnitude of the path update
  const double fA;         ///< Time step in the time discretization
  std::vector<double>
      fPath;  ///< Vector containing the current path configuration: it is initialized to zero
  std::vector<std::vector<double>> fResult;  ///< Matrix containing the sampled configurations
  std::vector<std::vector<std::vector<double>>>
      fBootStrapSet;  ///< Set of matrices containing the bootstrap copies of the initial fResult
                      ///< set

  /// Method to plot the status during execution
  ///
  /// \param index index of the current step inside the loop
  /// \param N_max maximum index, corresponding to the loop ending
  void PrintStatus(int index, int N_max) const;

 public:
  Metropolis() = delete;

  /// Metropolis constructor
  ///
  /// Constructor of the Metropolis class
  /// \param N number of grid points for the discretized trajectory
  /// \param N_corr number of correlated configurations to skip before next acquisition
  /// \param N_cf number of configurations to sample
  /// \param N_bootstraps number of statistical bootstraps to perform when BootStrap function is
  /// called \param epsilon typical magnitude of the path update \param a time step in the time
  /// discretization
  Metropolis(int N, int N_corr, int N_cf, int N_bootstraps, double epsilon, double a);

  /// Metropolis destructor
  ~Metropolis();

  /// \return N, number of grid points for the discretized trajectory
  int GetN() const;

  /// \return Ncorr, number of correlated configurations to skip before next acquisition
  int GetNcorr() const;

  /// \return Ncf, number of configurations to sample
  int GetNcf() const;

  /// \return epsilon, typical magnitude of the path update
  double GetEpsilon() const;

  /// \return a, time step in the time discretization
  double GetA() const;

  /// \return current path configuration, as stored in the private attribute fPath
  std::vector<double> GetCurrentPath() const;

  /// \return set of sampled Montecarlo path configurations
  std::vector<std::vector<double>> GetCurrentResult() const;

  /// Print the current path configuration on standard output
  void PrintCurrentPath() const;

  /// Action
  ///
  /// Definition of the part of the action of the physical system which
  /// involves the point j of the grid. Implementation is in the file
  /// CUSTOM.h, which allows for user modification.
  /// Boundary conditions must be explicitely or implicitely specified for
  /// the start-end points corresponding to index=0 and index=N-1, N being
  /// the number of grid points in the discretized path.
  /// \param j index of the time position in the discretized path
  double S(int j) const;

  /// Update the current path once
  ///
  /// Implementation of the Metropolis algorithm for a single sweep and
  /// update over the 1D lattice.
  /// \return acceptance ratio of the updates
  double UpdateCurrentPath();

  /// Run the Metropolis algorithm
  ///
  /// Perform a complete run of the Metropolis algorithm and obtain the
  /// configuration sample.
  void RunMetropolis();

  /// Statistical bootstrap
  ///
  /// Perform a statistical bootstrap to obtain the NBootStraps copies of the
  /// sampled configurations.
  void BootStrap();

  /// Evaluate observable gamma
  ///
  /// Define and evaluate the observable gamma whose expectation value is computed in
  /// Metropolis::ComputeEstimators. Implementation is in the file
  /// CUSTOM.h, which allows for user modification.
  /// \param path path configuration over which the observable is defined
  /// \param n if the observable depends on time, n is the index corresponding to
  ///          time t=na, where a is the time step.
  double EvaluateGamma(const std::vector<double>& path, int n) const;

  /// Evaluate propagator
  ///
  /// Define and evaluate the propagator whose expectation value is computed in
  /// Metropolis::ComputeEnergyEstimators.
  /// \param path path configuration over which the propagator is computed
  /// \param n integer index corresponding to time t=na, where a is the time step.
  double EvaluatePropagator(const std::vector<double>& path, int n) const;

  /// Compute gamma estimators
  ///
  /// Compute and print Montecarlo estimators (expectation value and 1-sigma
  /// uncertainty) for the observable defined in Metropolis::EvaluateGamma.
  void ComputeGammaEstimators();

  /// Compute energy estimators
  ///
  /// Compute and print on file the Montecarlo estimators for the energy gap
  /// between ground state and first excitated state. Uncertainty is estimated
  /// by building a set of bootstrap copies of the Montecarlo ensemble and evaluating
  /// the associated standard deviation.
  /// \param filename name of the output file
  void ComputeEnergyEstimators(std::string filename);

  /// Compute binned energy estimators
  ///
  /// Compute and print on file the Montecarlo estimators for the energy gap
  /// between ground state and first excitated state. Uncertainty is estimated
  /// by building a set of bootstrap copies of the Montecarlo ensemble and evaluating
  /// the associated standard deviation. In this method, before the bootstrap ensemble
  /// is built, a binning procedure is performed on the set of Montecarlo configurations.
  /// \param bin_width integer bin width, i.e. number of configurations to average to build a bin
  /// \param output_name name of the output file
  void ComputeBinnedEnergyEstimators(int bin_width, std::string output_name) const;
};

#endif