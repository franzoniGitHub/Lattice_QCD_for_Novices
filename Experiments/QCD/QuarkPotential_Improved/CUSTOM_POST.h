////////////////////////////////////////////////////////////////////////
/// \file CUSTOM_POST.h
/// \brief Header file to define a custom analysis function
///
/// Define here the function of the link variables whose expectation value should be computed.
/// Please, go through all the comments in the code, as you may find them helpful.
/// \see Type
///
////////////////////////////////////////////////////////////////////////
#ifndef CUSTOM_POST_H
#define CUSTOM_POST_H

// Write here the function of the link variables whose expectation value should be computed
// Please, go through all the comments here, as you may find them helpful

void Metropolis::ComputeCustom() const
{
  //############## IMPORTANT #################
  // Set the multiplicity parameter, which indicates how many equivalent terms are summed over in
  // each space-time point: if just one, set to 1.
  double multiplicity = 1.;
  //##########################################

  // Just some output and preliminary initializations
  std::cout << "Computing custom function expectation value..\nProgress %: " << std::flush;
  std::vector<int> n = fPath.GetNCells();  // n is now the vector containing the lattice dimensions
                                           // in the 4 directions

  double estimator = 0.;  // Estimator for the expectation value of the function
  double error = 0.;      // Estimator for the error, computed here as a standard deviation
  double square_estimator =
      0.;  // Estimator for the expectation value of the square of the function
  double old_estimator = 0.;

  for (int i = 0; i < fNcf; i++) {  // Loop over the lattice configurations: there are fNcf of them
    PrintStatus(i, fNcf);
    old_estimator = estimator;
    Path U = fResult[i];  // U is the whole lattice configuration in the i-th configuration
    // Here below, sweep over the lattice: modify from here, as you may or may not need a swipe over
    // the lattice
    for (int i0 = 0; i0 < n[0]; i0++) {
      for (int i1 = 0; i1 < n[1]; i1++) {
        for (int i2 = 0; i2 < n[2]; i2++) {
          for (int i3 = 0; i3 < n[3]; i3++) {
            for (int mu = 0; mu < 4; mu++) {     // this loop is over the polarizations mu of U_mu
              my4Vector x({i0, i1, i2, i3}, n);  // This defines the space-time point x
                                                 // To access to U_mu(x), just use U(x,mu)
              estimator += 0.;  // Modify this function
            }
          }
        }
      }
    }
    square_estimator += std::pow(
        (estimator - old_estimator) / (double)(n[0] * n[1] * n[2] * n[3] * multiplicity), 2.0);
  }
  cout << "\nEnd of statistics computation\n";
  estimator /= (double)(n[0] * n[1] * n[2] * n[3] * multiplicity * fNcf);
  square_estimator /= (double)fNcf;
  error = std::sqrt((square_estimator - std::pow(estimator, 2.0)) / fNcf);

  std::cout << "Results:\n<F[U_x_mu]> =   " << estimator << "  +/-  " << error << std::endl;
}

#endif
