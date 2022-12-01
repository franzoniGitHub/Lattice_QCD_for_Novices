#include <algorithm>
#include <armadillo>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <vector>
#include "my4Vector.h"
#include "Path.h"
#include "Metropolis.h"

using namespace arma;

/************************ Private Methods ***************************/

// Auxiliary method to compute action terms which don't depend on the updated
// link: implemented for the standard Wilson action
cx_dmat Metropolis::Gamma(const my4Vector& x, int mu) const
{
  cx_dmat result(3, 3, fill::zeros);
  for (int nu = 0; nu < 4; nu++) {
    if (nu != mu) {
      result = result +
               fPath(x.Offset(1, mu), nu) * fPath(x.Offset(1, nu), mu).t() * fPath(x, nu).t() +
               fPath(x.Offset(1, mu).Offset(-1, nu), nu).t() * fPath(x.Offset(-1, nu), mu).t() *
                   fPath(x.Offset(-1, nu), nu);
    }
  }
  return result;
}

// Auxiliary method to compute action terms which don't depend on the updated
// link: implemented for improved Wilson action
cx_dmat Metropolis::GammaImproved(const my4Vector& x, int mu) const
{
  cx_dmat result(3, 3, fill::zeros);
  for (int nu = 0; nu < 4; nu++) {
    if (nu != mu) {
      result = result +
               fPath(x.Offset(1, mu), mu) * fPath(x.Offset(2, mu), nu) *
                   fPath(x.Offset(1, mu).Offset(1, nu), mu).t() * fPath(x.Offset(1, nu), mu).t() *
                   fPath(x, nu).t() +
               fPath(x.Offset(1, mu), mu) * fPath(x.Offset(2, mu).Offset(-1, nu), nu).t() *
                   fPath(x.Offset(1, mu).Offset(-1, nu), mu).t() * fPath(x.Offset(-1, nu), mu).t() *
                   fPath(x.Offset(-1, nu), nu) +
               fPath(x.Offset(1, mu), nu) * fPath(x.Offset(1, mu).Offset(1, nu), nu) *
                   fPath(x.Offset(2, nu), mu).t() * fPath(x.Offset(1, nu), nu).t() *
                   fPath(x, nu).t() +
               fPath(x.Offset(1, mu).Offset(-1, nu), nu).t() *
                   fPath(x.Offset(1, mu).Offset(-2, nu), nu).t() * fPath(x.Offset(-2, nu), mu).t() *
                   fPath(x.Offset(-2, nu), nu) * fPath(x.Offset(-1, nu), nu) +
               fPath(x.Offset(1, mu), nu) * fPath(x.Offset(1, nu), mu).t() *
                   fPath(x.Offset(-1, mu).Offset(1, nu), mu).t() * fPath(x.Offset(-1, mu), nu).t() *
                   fPath(x.Offset(-1, mu), mu) +
               fPath(x.Offset(1, mu).Offset(-1, nu), nu).t() * fPath(x.Offset(-1, nu), mu).t() *
                   fPath(x.Offset(-1, mu).Offset(-1, nu), mu).t() *
                   fPath(x.Offset(-1, mu).Offset(-1, nu), nu) * fPath(x.Offset(-1, mu), mu);
    }
  }
  return result;
}

// S computes the action associated with the link U_mu(x)
// WARNING: gamma MUST be the result of Gamma(x, mu)
double Metropolis::S(const my4Vector& x,
                     int mu,
                     const cx_dmat& gamma_x_mu,
                     const cx_dmat& gamma_improved_x_mu) const
{
  if (fImproved)
    return (-fBetaTilde / 3.) *
           ((5. / (3. * std::pow(fU0, 4.))) * trace(fPath(x, mu) * gamma_x_mu).real() -
            (1. / (12. * std::pow(fU0, 6.))) * trace(fPath(x, mu) * gamma_improved_x_mu).real());
  else
    return (-fBeta / 3.) * trace(fPath(x, mu) * gamma_x_mu).real();
}

// Auxiliary method to print a status bar when performing time demanding loops
void Metropolis::PrintStatus(int index, int n_max) const
{
  double percentage = (((double)index) * 100.0) / (double)n_max;
  if (fmod(percentage, 5.) < 1.e-5) std::cout << percentage << " - " << std::flush;
  if (index == n_max - 1) std::cout << 100 << std::endl << std::flush;
}

// Auxiliary method to print the experiment settings as a header of an output file
void Metropolis::PrintSettingsOnFile(std::ofstream& file) const
{
  std::vector<int> n = fPath.GetNCells();
  file << "########## SETTINGS RECAP ##########\n";
  file << "Lattice dimension: " << n[0] << "  " << n[1] << "  " << n[2] << "  " << n[3]
       << std::endl;
  file << "Number of SU3: " << fNofSU3 << std::endl;
  file << "Number of correlated configurations to skip: " << fNcorr << std::endl;
  file << "Number of updates on each link variable: " << fInnerCycles << std::endl;
  file << "Number of sampled configurations: " << fNcf << std::endl;
  file << "Grid spacing: " << fA << std::endl;
  file << "Beta: " << fBeta << std::endl;
  file << "Beta_tilde: " << fBetaTilde << std::endl;
  file << "u0 coefficient: " << fU0 << std::endl;
  file << "Improved? " << std::boolalpha << fImproved << std::noboolalpha << std::endl;
  file << "####################################\n\n";
}

/************************ Public Methods ***************************/

// Constructor from the SETTINGS_EXP.h parameters
Metropolis::Metropolis(std::vector<int> N,
                       std::vector<int> integer_params,
                       std::vector<double> floating_params,
                       bool isimproved)
    : fNofSU3(integer_params[0])
    , fNcorr(integer_params[1])
    , fInnerCycles(integer_params[2])
    , fNcf(integer_params[3])
    , fA(floating_params[0])
    , fBeta(floating_params[1])
    , fBetaTilde(floating_params[2])
    , fU0(floating_params[3])
    , fEpsilon(floating_params[4])
    , fImproved(isimproved)
    , fPath(N)
{
  if ((integer_params.size() != 4) || (floating_params.size() != 5)) throw 1;
  for (int i = 0; i < fNcf; i++) fResult.push_back(fPath);
  for (int i = 0; i < 2 * fNofSU3; i++) fSetOfSU3.push_back(cx_dmat(3, 3, fill::zeros));
}

// Constructor from an input file
Metropolis::Metropolis(std::string infile) : fPath()
{
  std::ifstream file_input(infile);
  if (!file_input) {
    cout << "ERROR while opening the input file.\n";
    throw 1;
  }
  file_input >> fNofSU3 >> fNcorr >> fInnerCycles >> fNcf;
  file_input >> fA >> fBeta >> fBetaTilde >> fU0 >> fEpsilon;
  file_input >> fImproved;
  std::vector<int> n = {0, 0, 0, 0};
  file_input >> n[0] >> n[1] >> n[2] >> n[3];
  fPath.Reshape(n);
  for (int i = 0; i < fNcf; i++) fResult.push_back(fPath);
  double real = 0., imag = 0.;
  for (int index = 0; index < fNcf; index++) {
    for (int i0 = 0; i0 < n[0]; i0++) {
      for (int i1 = 0; i1 < n[1]; i1++) {
        for (int i2 = 0; i2 < n[2]; i2++) {
          for (int i3 = 0; i3 < n[3]; i3++) {
            for (int mu = 0; mu < 4; mu++) {
              for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                  my4Vector x({i0, i1, i2, i3}, n);
                  file_input >> real >> imag;
                  if (file_input) {
                    std::complex<double> matrix_element(real, imag);
                    my4Vector x({i0, i1, i2, i3}, n);
                    fResult[index].Set(x, mu)(i, j) = matrix_element;
                  } else {
                    std::cout << "ERROR while reading the lattice configurations from file.\n"
                              << std::flush;
                    throw 1;
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  file_input.close();
}

// Destructor
Metropolis::~Metropolis()
{
}

// Getters
std::vector<int> Metropolis::GetNCells() const
{
  return fPath.GetNCells();
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
double Metropolis::GetBeta() const
{
  return fBeta;
}
double Metropolis::GetU0() const
{
  return fU0;
}
std::vector<cx_dmat> Metropolis::GetSetOfSU3() const
{
  return fSetOfSU3;
}
Path Metropolis::GetCurrentPath() const
{
  return fPath;
}
std::vector<Path> Metropolis::GetCurrentResult() const
{
  return fResult;
}
bool Metropolis::IsImproved() const
{
  return fImproved;
}

// Print current path on standard output
void Metropolis::PrintPathOnScreen() const
{
  fPath.Print();
}

// Print current path on file
void Metropolis::PrintPathOnFile(std::string filename) const
{
  typedef std::numeric_limits<double> dbl;
  const auto default_precision = std::cout.precision();
  std::cout << std::scientific << std::setprecision(dbl::digits10);

  std::vector<int> n = fPath.GetNCells();
  std::ofstream file_path(filename);
  for (int i0 = 0; i0 < n[0]; i0++) {
    for (int i1 = 0; i1 < n[1]; i1++) {
      for (int i2 = 0; i2 < n[2]; i2++) {
        for (int i3 = 0; i3 < n[3]; i3++) {
          for (int mu = 0; mu < 4; mu++) {
            for (int i = 0; i < 3; i++) {
              for (int j = 0; j < 3; j++) {
                my4Vector x({i0, i1, i2, i3}, n);
                file_path << i0 << "  " << i1 << "  " << i2 << "  " << i3 << "  " << mu << "  " << i
                          << "  " << j << "  " << fPath(x, mu)(i, j).real() << "  "
                          << fPath(x, mu)(i, j).imag() << std::endl;
              }
            }
          }
        }
      }
    }
  }
  file_path.close();
  std::cout << std::defaultfloat << std::setprecision(default_precision);
}

// Clear result: bring the result vector to the default one
void Metropolis::ClearResult()
{
  for (int i = 0; i < fNcf; i++) fResult[i] = Path(fPath.GetNCells());
}

// Print settings and Montecarlo configurations on file:
// by default, the print mode is silent, which is the required option for later postprocessing
// (matrix elements are printed with maximum accuracy, no headers); as an option, 'V' may be
// specified for a verbose print, which is more user readable, yet less accurate (it saves disk
// space, though).
void Metropolis::PrintAllOnFile(std::string filename, char mode) const
{
  std::cout << "Printing the lattice configurations on file..\n" << std::flush;
  std::vector<int> n = fPath.GetNCells();
  std::ofstream file_result(filename);
  if (mode == 'V') {  // Verbose option
    PrintSettingsOnFile(file_result);
    file_result << "Sampled configurations will follow:\n";
    file_result << "ith-configuration; (x, y, x, t, mu) grid position and mu index, (l,m) matrix "
                   "component, a+ib complex matrix element\n";
    for (int index = 0; index < fNcf; index++) {
      for (int i0 = 0; i0 < n[0]; i0++) {
        for (int i1 = 0; i1 < n[1]; i1++) {
          for (int i2 = 0; i2 < n[2]; i2++) {
            for (int i3 = 0; i3 < n[3]; i3++) {
              for (int mu = 0; mu < 4; mu++) {
                for (int i = 0; i < 3; i++) {
                  for (int j = 0; j < 3; j++) {
                    my4Vector x({i0, i1, i2, i3}, n);
                    file_result << index << "; (" << i0 << ", " << i1 << ", " << i2 << ", " << i3
                                << ", " << mu << "), (" << i << ", " << j << "), "
                                << fResult[index](x, mu)(i, j).real() << "+i"
                                << fResult[index](x, mu)(i, j).imag() << std::endl;
                  }
                }
              }
            }
          }
        }
      }
    }
  } else {  // Default, silent option
    typedef std::numeric_limits<double> dbl;
    const auto default_precision = file_result.precision();
    file_result << std::scientific << std::setprecision(dbl::digits10);

    file_result << fNofSU3 << std::endl;
    file_result << fNcorr << std::endl;
    file_result << fInnerCycles << std::endl;
    file_result << fNcf << std::endl;
    file_result << fA << std::endl;
    file_result << fBeta << std::endl;
    file_result << fBetaTilde << std::endl;
    file_result << fU0 << std::endl;
    file_result << fEpsilon << std::endl;
    file_result << fImproved << std::endl;
    file_result << n[0] << "  " << n[1] << "  " << n[2] << "  " << n[3] << std::endl;
    for (int index = 0; index < fNcf; index++) {
      for (int i0 = 0; i0 < n[0]; i0++) {
        for (int i1 = 0; i1 < n[1]; i1++) {
          for (int i2 = 0; i2 < n[2]; i2++) {
            for (int i3 = 0; i3 < n[3]; i3++) {
              for (int mu = 0; mu < 4; mu++) {
                for (int i = 0; i < 3; i++) {
                  for (int j = 0; j < 3; j++) {
                    my4Vector x({i0, i1, i2, i3}, n);
                    file_result << fResult[index](x, mu)(i, j).real() << " "
                                << fResult[index](x, mu)(i, j).imag() << " ";
                  }
                }
              }
            }
          }
        }
      }
      file_result << std::endl;
    }
    file_result << std::defaultfloat << std::setprecision(default_precision);
  }
  file_result.close();
}

// Randomize to get the set of random SU3 matrices
void Metropolis::RandomizeSU3()
{
  std::random_device rd;
  std::mt19937_64 generator(rd());
  std::uniform_real_distribution<> uni(-1.0, 1.0);
  // Cycle over the set of matrices
  for (int i = 0; i < fNofSU3; i++) {
    fSetOfSU3[i].zeros();
    // Update the entries and build random matrices
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        fSetOfSU3[i](j, k) = std::complex<double>(uni(generator), uni(generator));
      }
    }
    // Make the previous matrices hermitian
    fSetOfSU3[i] = 0.5 * (fSetOfSU3[i] + fSetOfSU3[i].t());
    // Build SU3 by complex matrix exponentiation of the previous hermitian matrices
    // and then dividing by determinant^(1/3)
    fSetOfSU3[i] = expmat((fEpsilon * std::complex<double>(0., 1.)) * fSetOfSU3[i]);
    fSetOfSU3[i] = fSetOfSU3[i] / std::pow(det(fSetOfSU3[i]), 1. / 3.);
    // Append the adjoints
    fSetOfSU3[i + fNofSU3] = fSetOfSU3[i].t();
  }
}

// Compute plaquette and rectangle expectation values
void Metropolis::ComputePlaquetteRectangle() const
{
  std::cout << "Computing 1x1 and 2x1 plaquettes expectation values..\nProgress %: " << std::flush;

  std::vector<int> n = fPath.GetNCells();
  std::vector<double> estimators = {0., 0.};  // First= axa plaquette; Second=ax2a rectangle
  std::vector<double> errors = {0., 0.};      // Same as above
  std::vector<double> square_estimators = {0., 0.};
  std::vector<double> old_estimators = {0., 0.};
  for (int i = 0; i < fNcf; i++) {
    PrintStatus(i, fNcf);
    old_estimators[0] = estimators[0];
    old_estimators[1] = estimators[1];
    for (int i0 = 0; i0 < n[0]; i0++) {
      for (int i1 = 0; i1 < n[1]; i1++) {
        for (int i2 = 0; i2 < n[2]; i2++) {
          for (int i3 = 0; i3 < n[3]; i3++) {
            for (int mu = 0; mu < 4; mu++) {
              for (int nu = 0; nu < mu; nu++) {
                my4Vector x({i0, i1, i2, i3}, n);
                estimators[0] += WilsonLoop(1, 1, mu, nu, x, i);
                estimators[1] += WilsonLoop(2, 1, mu, nu, x, i);
              }
            }
          }
        }
      }
    }
    square_estimators[0] += std::pow(
        (estimators[0] - old_estimators[0]) / (double)(n[0] * n[1] * n[2] * n[3] * 6.), 2.0);
    square_estimators[1] += std::pow(
        (estimators[1] - old_estimators[1]) / (double)(n[0] * n[1] * n[2] * n[3] * 6.), 2.0);
  }
  std::cout << "End of statistics computation\n";
  estimators[0] /= (double)(n[0] * n[1] * n[2] * n[3] * 6. * fNcf);
  estimators[1] /= (double)(n[0] * n[1] * n[2] * n[3] * 6. * fNcf);

  square_estimators[0] /= (double)fNcf;
  square_estimators[1] /= (double)fNcf;

  errors[0] = std::sqrt((square_estimators[0] - pow(estimators[0], 2.0)) / fNcf);
  errors[1] = std::sqrt((square_estimators[1] - pow(estimators[1], 2.0)) / fNcf);

  std::cout << "Results:\nsquare plaquette =   " << estimators[0] << "  +/-  " << errors[0]
            << std::endl;
  std::cout << "rectangular plaquette =   " << estimators[1] << "  +/-  " << errors[1] << std::endl;
}

// Compute RxT Wilson Loops
void Metropolis::ComputeRxTWilsonLoops() const
{
  std::cout << "Computing RxT Wilson loops..\nProgress %: " << std::flush;
  std::vector<int> n = fPath.GetNCells();
  // Define the max loop dimensions in the R and T directions as (approximately) half the grid dimensions
  const int nR = (int)(std::min({n[0], n[1], n[2]}) / 2.);
  const int nT = (int)(n[3] / 2.);
  // Define the matrices containing the results for the RxT loop dimension
  std::vector<std::vector<double>> estimators, errors, square_estimators, old_estimators;
  std::vector<double> inner;
  for (int R = 0; R <= nR; R++) inner.push_back(0.);
  for (int T = 0; T <= nT; T++) estimators.push_back(inner);
  errors = estimators;
  square_estimators = estimators;
  old_estimators = estimators;

  // Loop over the configurations and the lattice to get the Montecarlo estimators for the loops
  for (int i = 0; i < fNcf; i++) {
    PrintStatus(i, fNcf);
    old_estimators = estimators;
    for (int i0 = 0; i0 < n[0]; i0++) {
      for (int i1 = 0; i1 < n[1]; i1++) {
        for (int i2 = 0; i2 < n[2]; i2++) {
          for (int i3 = 0; i3 < n[3]; i3++) {
            my4Vector x({i0, i1, i2, i3}, n);
            for (int T = 1; T <= nT; T++) {
              for (int R = 1; R <= nR; R++) {
                for (int space_dir = 0; space_dir < 3; space_dir++) {
                  estimators[T][R] +=
                      WilsonLoop(T, R, 3, space_dir, x, i) + WilsonLoop(R, T, space_dir, 3, x, i);
                }
              }
            }
          }
        }
      }
    }
    for (int T = 1; T <= nT; T++) {
      for (int R = 1; R <= nR; R++) {
        square_estimators[T][R] += std::pow(
            (estimators[T][R] - old_estimators[T][R]) / (double)(n[0] * n[1] * n[2] * n[3] * 6.),
            2.0);
      }
    }
  }

  // File output: first print the loop estimators
  std::ofstream file_loop_output("RXT_loops_file.dat");
  PrintSettingsOnFile(file_loop_output);
  file_loop_output << "rxt planar Wilson loop expectation values (the format is -> value:error)\n";
  file_loop_output << "Vertical: t/a      Horizontal: r/a\n";
  for (int R = 1; R <= nR; R++) file_loop_output << "\t" << R;
  file_loop_output << std::endl;
  for (int T = 1; T <= nT; T++) {
    file_loop_output << T;
    for (int R = 1; R <= nR; R++) {
      estimators[T][R] /= n[0] * n[1] * n[2] * n[3] * 6. * fNcf;
      errors[T][R] = sqrt((square_estimators[T][R] / fNcf - pow(estimators[T][R], 2.0)) / fNcf);
      file_loop_output << "\t" << estimators[T][R] << ":" << errors[T][R];
    }
    file_loop_output << std::endl;
  }
  file_loop_output.close();
  std::cout << "Printing on file \"RXT_loops_file.dat\"\n";

  // File output: second print the W(r,t)/W(r,t+a) estimators
  std::ofstream file_potential_output("RXT_potential_file.dat");
  PrintSettingsOnFile(file_potential_output);
  file_potential_output << "Ratio W(r,t)/W(r,t+a) (the format is -> value:error)\n";
  file_potential_output << "Vertical: t/a        Horizontal: r/a\n";
  for (int R = 1; R <= nR; R++) file_potential_output << "\t" << R;
  file_potential_output << std::endl;
  for (int T = 1; T < nT; T++) {
    file_potential_output << T;
    for (int R = 1; R <= nR; R++) {
      file_potential_output << "\t" << estimators[T][R] / estimators[T + 1][R] << ":"
                            << std::pow(std::pow(errors[T][R] / estimators[T][R], 2.) +
                                            std::pow(errors[T + 1][R] / estimators[T + 1][R], 2.),
                                        0.5) *
                                   std::abs(estimators[T][R] / estimators[T + 1][R]);
    }
    file_potential_output << std::endl;
  }
  std::cout << "Printing on file \"RXT_potential_file.dat\"\n";
  file_potential_output.close();

  // File output: third print the asymptotic W(r,t)/W(r,t+a) estimators, i.e. the potential
  // estimators
  //              and propagate the uncertainty.
  std::ofstream file_potential_plot("RXT_potential_plot_file.dat");
  PrintSettingsOnFile(file_potential_plot);
  file_potential_plot << "Potential aV(r)\n";
  file_potential_plot << "First column: r/a    Second column: aV(r)   Third column: sigma(aV(r))\n";
  for (int R = 1; R <= nR; R++) {
    file_potential_plot << R << "\t" << estimators[nT - 1][R] / estimators[nT][R] << "\t"
                        << std::pow(std::pow(errors[nT - 1][R] / estimators[nT - 1][R], 2.) +
                                        std::pow(errors[nT][R] / estimators[nT][R], 2.),
                                    0.5) *
                               std::abs(estimators[nT - 1][R] / estimators[nT][R])
                        << std::endl;
  }
  file_potential_plot.close();
  std::cout << "Printing on file \"RXT_potential_plot_file.dat\"\n";
}

#include "../CUSTOM_POST.h"  //This defines the custom statistics function

// Compute and return gauge-covariant derivative summed over all directions on the i-th result
// element
Path Metropolis::GaugeDerivative(int i) const
{
  Path outPath(fPath.GetNCells());
  Path U = fResult[i];
  std::vector<int> n = fPath.GetNCells();
  for (int i0 = 0; i0 < n[0]; i0++) {
    for (int i1 = 0; i1 < n[1]; i1++) {
      for (int i2 = 0; i2 < n[2]; i2++) {
        for (int i3 = 0; i3 < n[3]; i3++) {
          for (int mu = 0; mu < 4; mu++) {
            cx_dmat sum(3, 3, fill::zeros);
            my4Vector x({i0, i1, i2, i3}, n);
            for (int rho = 0; rho < 4; rho++) {
              sum += (1. / std::pow(fU0 * fA, 2.)) *
                     (U(x, rho) * U(x.Offset(1, rho), mu) * U(x.Offset(1, mu), rho).t() -
                      2. * fU0 * fU0 * U(x, mu) +
                      U(x.Offset(-1, rho), rho).t() * U(x.Offset(-1, rho), mu) *
                          U(x.Offset(-1, rho).Offset(1, mu), rho));
            }
            outPath.Set(x, mu) = sum;
          }
        }
      }
    }
  }
  return outPath;
}

// Smear fResult NTimes with a smearing parameter smearing_par
void Metropolis::SpatialSmearing(int Ntimes, double smearing_par)
{
  std::vector<int> n = fPath.GetNCells();
  std::cout << "Spatial smearing of the link variables in the results (" << Ntimes
            << " times)..\nProgress %: ";
  int iteration = 0;
  for (int i_smear = 0; i_smear < Ntimes; i_smear++) {
    for (int i = 0; i < fNcf; i++) {
      Path gauge_der = GaugeDerivative(i);
      for (int i0 = 0; i0 < n[0]; i0++) {
        for (int i1 = 0; i1 < n[1]; i1++) {
          for (int i2 = 0; i2 < n[2]; i2++) {
            for (int i3 = 0; i3 < n[3]; i3++) {
              for (int mu = 0; mu < 3; mu++) {
                my4Vector x({i0, i1, i2, i3}, n);
                fResult[i].Set(x, mu) =
                    fResult[i](x, mu) + smearing_par * fA * fA * gauge_der(x, mu);
              }
            }
          }
        }
      }
      PrintStatus(iteration, fNcf * Ntimes);
      iteration++;
    }
  }
}

// Update the current path: it returns the acceptance ratio for the update
double Metropolis::UpdateCurrentPath()
{
  // Define the pseudo-random Mersenne-Twister generators
  std::random_device rd;
  std::mt19937_64 generator1(rd());
  std::mt19937_64 generator2(rd());
  // Define the probability distributions
  std::uniform_real_distribution<> uni_01(0.0, 1.0);
  std::uniform_int_distribution<> uni_int(0, 2 * fNofSU3 - 1);
  double accepted = 0.0;
  std::vector<int> n = fPath.GetNCells();
  // Seep over the lattice and do the update
  for (int i0 = 0; i0 < n[0]; i0++) {
    for (int i1 = 0; i1 < n[1]; i1++) {
      for (int i2 = 0; i2 < n[2]; i2++) {
        for (int i3 = 0; i3 < n[3]; i3++) {
          for (int mu = 0; mu < 4; mu++) {
            my4Vector x({i0, i1, i2, i3}, n);
            // Compute the action independent on the link at (x,mu)
            cx_dmat gamma_x_mu = Gamma(x, mu);
            cx_dmat gamma_improved_x_mu(3, 3, fill::zeros);
            if (fImproved) gamma_improved_x_mu = GammaImproved(x, mu);
            // Do fInnerCycles updates before going to the next site
            for (int inner = 0; inner < fInnerCycles; inner++) {
              cx_dmat old_link_x_mu = fPath(x, mu);
              double old_S_x_mu = S(x, mu, gamma_x_mu, gamma_improved_x_mu);
              int index = uni_int(generator1);
              fPath.Set(x, mu) = fSetOfSU3[index] * fPath(x, mu);
              double deltaS = S(x, mu, gamma_x_mu, gamma_improved_x_mu) - old_S_x_mu;
              // Accept or reject the update, depending on the sign of deltaS
              if ((deltaS < 0) || (std::exp(-deltaS) > uni_01(generator2)))
                accepted += 1.0;
              else
                fPath.Set(x, mu) = old_link_x_mu;
            }
          }
        }
      }
    }
  }
  return accepted / (double)(n[0] * n[1] * n[2] * n[3] * 4 * fInnerCycles);
}

// Run the Metropolis algorithm on the current path
void Metropolis::RunMetropolis()
{
  std::cout << "Randomization of the SU3 matrices...\n";
  RandomizeSU3();
  std::cout << "Metropolis is running..\nGrid thermalization...\nProgress %: " << std::flush;
  double avg = 0.0;
  int counter = 0;
  fResult.clear();

  for (int i = 0; i < 10 * fNcorr; i++) {  // thermalize the path
    PrintStatus(i, 10 * fNcorr);
    UpdateCurrentPath();
  }
  std::cout << "Generating configurations...\nProgress %: " << std::flush;
  for (int i = 0; i < fNcf; i++) {  // repeat the following Ncf times
    PrintStatus(i, fNcf);
    fResult.push_back(fPath);           // save current path
    for (int j = 0; j < fNcorr; j++) {  // discard fNcorr paths before saving again
      avg += UpdateCurrentPath();       // compute the average acceptance ratio
      counter++;
    }
  }

  std::cout << "Metropolis has finished. The avg acceptance level is: " << avg / counter
            << std::endl;
}

// WILSON LOOP: to compute N_mu x N_nu Wilson loops around position x in the mu-nu plane using the
// j-th path configuration
double Metropolis::WilsonLoop(int N_mu, int N_nu, int mu, int nu, const my4Vector& x, int j) const
{
  // lower horizontal segment in the mu direction
  cx_dmat LowerMu(3, 3, fill::eye);
  for (int i = 0; i < N_mu; i++) LowerMu = LowerMu * fResult[j](x.Offset(i, mu), mu);
  // right vertical segment in the nu direction
  cx_dmat RightNu(3, 3, fill::eye);
  for (int i = 0; i < N_nu; i++)
    RightNu = RightNu * fResult[j](x.Offset(N_mu, mu).Offset(i, nu), nu);
  // upper horizontal segment "backward" in the mu direction
  cx_dmat UpperMu(3, 3, fill::eye);
  for (int i = 1; i <= N_mu; i++)
    UpperMu = UpperMu * (fResult[j](x.Offset(N_nu, nu).Offset(N_mu - i, mu), mu).t());
  // left vertical segment "downward" in the nu direction
  cx_dmat LeftNu(3, 3, fill::eye);
  for (int i = 1; i <= N_nu; i++) LeftNu = LeftNu * (fResult[j](x.Offset(N_nu - i, nu), nu).t());
  // combine the segments into a Wilson loop
  return (1. / 3.) * ((trace(LowerMu * RightNu * UpperMu * LeftNu)).real());
}

// Method to select the analysis to compute using the enum class Type
void Metropolis::ComputeStatistics(Type type) const
{
  switch (type) {
    case Type::PlaquetteRectangle: {
      ComputePlaquetteRectangle();
      break;
    }
    case Type::QuarkPotential: {
      ComputeRxTWilsonLoops();
      break;
    }
    case Type::Custom: {
      ComputeCustom();
      break;
    }
    default:
      std::cout << "No statistics computed, wrong predefined type\n";
  }
}
