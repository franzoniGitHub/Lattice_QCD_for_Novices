#include <chrono>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>
#include "Metropolis.h"
#include "../CUSTOM.h"  // Include here the implementations of S, EvaluateGamma and EvaluatePropagator

/************************ Private Methods ***************************/

// Print status
void Metropolis::PrintStatus(int index, int n_max) const
{
  double percentage = (((double)index) * 100.0) / (double)n_max;
  if (fmod(percentage, 5.) < 1.e-5) std::cout << percentage << " - " << std::flush;
  if (index == n_max - 1) std::cout << 100 << std::endl << std::flush;
}

/************************ Public Methods ***************************/

// Constructor
Metropolis::Metropolis(int N, int N_corr, int N_cf, int N_bootstraps, double epsilon, double a)
    : fN(N), fNcorr(N_corr), fNcf(N_cf), fNBootStraps(N_bootstraps), fEpsilon(epsilon), fA(a)
{
  std::vector<double> zero_path;
  for (int i = 0; i < N; i++) fPath.push_back(0.0);
  for (int i = 0; i < N_cf; i++) fResult.push_back(fPath);
  for (int i = 0; i < N_bootstraps; i++) fBootStrapSet.push_back(fResult);
}

// Destructor
Metropolis::~Metropolis()
{
}

// Getters
int Metropolis::GetN() const
{
  return fN;
}

int Metropolis::GetNcorr() const
{
  return fNcorr;
}

int Metropolis::GetNcf() const
{
  return fNcf;
}

double Metropolis::GetEpsilon() const
{
  return fEpsilon;
}

double Metropolis::GetA() const
{
  return fA;
}

std::vector<double> Metropolis::GetCurrentPath() const
{
  return fPath;
}

std::vector<std::vector<double>> Metropolis::GetCurrentResult() const
{
  return fResult;
}

// Print current path
void Metropolis::PrintCurrentPath() const
{
  for (int i = 0; i < fN; i++) std::cout << fPath[i] << "  ";
  std::cout << std::endl;
}

// Update the current path: it returns the acceptance ratio for the update
double Metropolis::UpdateCurrentPath()
{
  std::random_device rd;
  std::mt19937_64 generator1(rd());
  std::mt19937_64 generator2(rd());
  std::uniform_real_distribution<> uni_01(0.0, 1.0);
  std::uniform_real_distribution<> uni_epsilon(-fEpsilon, fEpsilon);
  double accepted = 0.0;
  for (int i = 0; i < fN; i++) {
    double old_path_i = fPath[i];
    double old_S_i = S(i);
    fPath[i] += uni_epsilon(generator1);
    double deltaS = S(i) - old_S_i;
    if ((deltaS > 0.) && (std::exp(-deltaS) < uni_01(generator2)))
      fPath[i] = old_path_i;
    else
      accepted += 1.0;
  }
  return accepted / (double)fN;
}

// Run the Metropolis algorithm to sample configurations
void Metropolis::RunMetropolis()
{
  std::cout << "Metropolis is running\nStatus (%): " << std::flush;
  double avg = 0.0;
  int counter = 0;
  fResult.clear();
  for (int i = 0; i < 5 * fNcorr; i++) UpdateCurrentPath();  // thermalize the path

  for (int i = 0; i < fNcf; i++) {  // repeat the following Ncf times
    PrintStatus(i, fNcf);
    fResult.push_back(fPath);           // save current path
    for (int j = 0; j < fNcorr; j++) {  // discard fNcorr paths before saving again
      avg += UpdateCurrentPath();       // compute the average acceptance ratio
      counter++;
    }
  }
  std::cout << "Metropolis has finished. The avg acceptance level is: " << avg / (double)counter
            << std::endl;
}

// Bootstrap fResult fNBootStraps times
void Metropolis::BootStrap()
{
  std::random_device rd;
  std::mt19937_64 generator(rd());
  std::uniform_int_distribution<> uni(0, fNcf - 1);
  fBootStrapSet[0] = fResult;
  for (int i = 1; i < fNBootStraps; i++) {
    for (int j = 0; j < fNcf; j++) {
      fBootStrapSet[i][j] = fResult[uni(generator)];
    }
  }
}

// Get gamma estimators by averaging on the Ncf configurations
void Metropolis::ComputeGammaEstimators()
{
  std::cout << "Computing the Gamma estimators by averaging on the Ncf configurations\n";
  std::cout << "time\tGamma(t)\tdGamma\n" << std::flush;
  // Computing the Gamma estimators on the Montecarlo ensemble
  for (int j = 0; j < fN; j++) {
    double partial_sum_g = 0.0;
    double partial_sum_g_square = 0.0;
    for (int k = 0; k < fNcf; k++) {
      partial_sum_g += EvaluateGamma(fResult[k], j);
      partial_sum_g_square += EvaluateGamma(fResult[k], j) * EvaluateGamma(fResult[k], j);
    }
    double avg = partial_sum_g / (double)fNcf;
    double error = std::sqrt((partial_sum_g_square / (double)fNcf - avg * avg) / (double)fNcf);
    std::cout << j * fA << "\t" << avg << "\t" << error << std::endl;
  }
}

// Get energy estimators using the set of bootstrap copies
void Metropolis::ComputeEnergyEstimators(std::string output_name)
{
  BootStrap();  // Build the set of bootstrap copies of the Montecarlo ensemble
  std::vector<std::vector<double>>
      prop_single;  // Vector of propagator {avg, error} for the fN times
  std::vector<std::vector<std::vector<double>>>
      prop_bootstrapset;  // Set of above vectors for every bootstrap copy
  std::cout << "Computing the propagator estimators by averaging on the Ncf configurations\n";
  std::cout << "Status (%): " << std::flush;
  int status = 0;
  // Computing the propagator estimators for every bootstrap copy and storing the result in
  // prop_bootstrapset
  for (int i = 0; i < fNBootStraps; i++) {
    for (int j = 0; j < fN; j++) {
      double partial_sum_p = 0.0;
      double partial_sum_p_square = 0.0;
      PrintStatus(status, fN * fNBootStraps);
      status++;
      for (int k = 0; k < fNcf; k++) {
        partial_sum_p += EvaluatePropagator(fBootStrapSet[i][k], j);
        partial_sum_p_square +=
            EvaluatePropagator(fBootStrapSet[i][k], j) * EvaluatePropagator(fBootStrapSet[i][k], j);
      }
      double avg = partial_sum_p / (double)fNcf;
      double error = std::sqrt((partial_sum_p_square / (double)fNcf - avg * avg) / (double)fNcf);
      prop_single.push_back({avg, error});
    }
    prop_bootstrapset.push_back(prop_single);
    prop_single.clear();
  }
  std::ofstream outfile(output_name);
  std::cout << "Computing deltaE estimators by averaging on the NBootStraps configurations:\n";
  std::cout << "time\tPropagator(t)\tdPropagator\tdeltaE(t)\tddeltaE\n" << std::flush;
  outfile << "time\tPropagator(t)\tdPropagator\tdeltaE(t)\tddeltaE\n";
  // Computing the energy estimators and printing all on file and standard output
  for (int n = 0; n < fN - 1; n++) {
    double sum_dE_n = 0.;
    double sum_dE2_n = 0.;
    for (int i = 0; i < fNBootStraps; i++) {
      sum_dE_n +=
          std::log(std::abs(prop_bootstrapset[i][n][0] / prop_bootstrapset[i][n + 1][0])) / fA;
      sum_dE2_n += std::pow(
          std::log(std::abs(prop_bootstrapset[i][n][0] / prop_bootstrapset[i][n + 1][0])) / fA, 2.);
    }
    std::cout << n * fA << "\t" << prop_bootstrapset[0][n][0] << "\t" << prop_bootstrapset[0][n][1];
    std::cout << "\t" << sum_dE_n / (double)fNBootStraps << "\t";
    std::cout << std::sqrt(sum_dE2_n / (double)fNBootStraps -
                           std::pow(sum_dE_n / (double)fNBootStraps, 2.))
              << std::endl
              << std::flush;

    outfile << n * fA << "\t" << prop_bootstrapset[0][n][0] << "\t" << prop_bootstrapset[0][n][1];
    outfile << "\t" << sum_dE_n / (double)fNBootStraps << "\t";
    outfile << std::sqrt(sum_dE2_n / (double)fNBootStraps -
                         std::pow(sum_dE_n / (double)fNBootStraps, 2.))
            << std::endl;
  }
}

// Get energy estimators performing a binning procedure which groups bin_width configurations before
// building the bootstrap set
void Metropolis::ComputeBinnedEnergyEstimators(int bin_width, std::string output_name) const
{
  // Initializations
  std::vector<double> zero_vector;
  std::vector<std::vector<double>> binned_propagator_set;
  for (int i = 0; i < fN; i++) zero_vector.push_back(0.);
  for (int i = 0; i < fNcf; i += bin_width) binned_propagator_set.push_back(zero_vector);

  // First, compute a binned average
  int binnedSize = binned_propagator_set.size();
  for (int j = 0; j < fN; j++) {
    int binnedIndex = 0;
    for (int k = 0; k < fNcf; k += bin_width) {
      if (k + bin_width >= fNcf) break;
      for (int i = 0; i < bin_width; i++) {
        binned_propagator_set[binnedIndex][j] += EvaluatePropagator(fResult[k + i], j);
      }
      binned_propagator_set[binnedIndex][j] /= (double)bin_width;
      binnedIndex++;
    }
  }

  // Second, create the bootstrap set of copies
  std::vector<std::vector<std::vector<double>>> binned_bootstrap_set;
  for (int i = 0; i < fNBootStraps; i++) binned_bootstrap_set.push_back(binned_propagator_set);
  std::random_device rd;
  std::mt19937_64 generator(rd());
  std::uniform_int_distribution<> uni(0, binnedSize - 1);
  for (int i = 1; i < fNBootStraps; i++) {
    for (int j = 0; j < binnedSize; j++) {
      binned_bootstrap_set[i][j] = binned_propagator_set[uni(generator)];
    }
  }

  std::vector<std::vector<double>>
      prop_single;  // Vector of propagator values {avg, error} for the fN times
  std::vector<std::vector<std::vector<double>>>
      prop_bootstrapset;  // Set of above vectors for every bootstrap copy
  std::cout << "Computing the propagator estimators by averaging on the binned configurations\n";
  std::cout << "Status (%): " << std::flush;
  int status = 0;

  // Third, compute the propagator estimators for every bootstrap copy and store the result in
  // prop_bootstrapset
  for (int i = 0; i < fNBootStraps; i++) {
    for (int j = 0; j < fN; j++) {
      double partial_sum_p = 0.0;
      double partial_sum_p_square = 0.0;
      PrintStatus(status, fN * fNBootStraps);
      status++;
      for (int k = 0; k < binnedSize; k++) {
        partial_sum_p += binned_bootstrap_set[i][k][j];
        partial_sum_p_square += binned_bootstrap_set[i][k][j] * binned_bootstrap_set[i][k][j];
      }
      double avg = partial_sum_p / (double)fNcf;
      double error = std::sqrt((partial_sum_p_square / (double)fNcf - avg * avg) / (double)fNcf);
      prop_single.push_back({avg, error});
    }
    prop_bootstrapset.push_back(prop_single);
    prop_single.clear();
  }

  // Finally, compute averages and errors on the bootstrap set
  std::ofstream outfile(output_name);
  std::cout << "Computing deltaE estimators by averaging on the NBootStraps configurations:\n";
  std::cout << "time\tPropagator(t)\tdPropagator\tdeltaE(t)\tddeltaE\n" << std::flush;
  outfile << "time\tPropagator(t)\tdPropagator\tdeltaE(t)\tddeltaE\n";

  for (int n = 0; n < fN - 1; n++) {
    double sum_dE_n = 0.;
    double sum_dE2_n = 0.;
    for (int i = 0; i < fNBootStraps; i++) {
      sum_dE_n +=
          std::log(std::abs(prop_bootstrapset[i][n][0] / prop_bootstrapset[i][n + 1][0])) / fA;
      sum_dE2_n += std::pow(
          std::log(std::abs(prop_bootstrapset[i][n][0] / prop_bootstrapset[i][n + 1][0])) / fA, 2.);
    }
    std::cout << n * fA << "\t" << prop_bootstrapset[0][n][0] << "\t" << prop_bootstrapset[0][n][1];
    std::cout << "\t" << sum_dE_n / (double)fNBootStraps << "\t";
    std::cout << std::sqrt(sum_dE2_n / (double)fNBootStraps -
                           std::pow(sum_dE_n / (double)fNBootStraps, 2.))
              << std::endl
              << std::flush;

    outfile << n * fA << "\t" << prop_bootstrapset[0][n][0] << "\t" << prop_bootstrapset[0][n][1];
    outfile << "\t" << sum_dE_n / (double)fNBootStraps << "\t";
    outfile << std::sqrt(sum_dE2_n / (double)fNBootStraps -
                         std::pow(sum_dE_n / (double)fNBootStraps, 2.))
            << std::endl;
  }
}
