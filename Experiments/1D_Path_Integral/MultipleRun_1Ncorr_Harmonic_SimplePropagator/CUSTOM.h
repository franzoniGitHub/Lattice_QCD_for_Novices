////////////////////////////////////////////////////////////////////////
/// \file CUSTOM.h
/// \brief Header file to define the action of the system and the observables
///        to be computed.
///
/// Define here the implementation of the action in Metropolis::S, that
/// of the observable gamma in Metropolis::EvaluateGamma and the one of
/// the propagator in Metropolis::EvaluatePropagator. Please, follow
/// the comments in the file for more specific indications.
///
////////////////////////////////////////////////////////////////////////

#ifndef CUSTOM_H
#define CUSTOM_H

// The method S defines the terms of the discretized action (imaginary time)
// involving the i-th site of the current path (the index i takes integer values on
// the range from 0 to fN). Please, remember to provide
// boundary conditions for the cases i-1<0 and i+1>=fN.\n
// NOTE: fPath is a vector containing the path configuration
double Metropolis::S(int i) const
{
  int i_next = (i + 1) % fN;  // next site position
  int i_previous = 0;         // previous site position
  if ((i - 1) < 0)
    i_previous = fN - 1;
  else
    i_previous = i - 1;
  return fA * fPath[i] * fPath[i] / 2.0 +
         fPath[i] * (fPath[i] - fPath[i_next] - fPath[i_previous]) / fA;
}

// Definition of the observable Gamma, whose expectation value is computed in
// the Metropolis::ComputeGammaEstimators.
// path: path on which the observable is evaluated
// n: time index taking values in the range of integers from 0 to fN -1.
double Metropolis::EvaluateGamma(const std::vector<double>& path, int n) const
{
  double partial_sum = 0.0;
  if (path.size() != fN) {
    std::cout << "ERROR: wrong path size in EvaluateGamma input\n";
    throw 1;
  }
  for (int k = 0; k < fN; k++) partial_sum += path[k] * path[(k + n) % fN];
  return partial_sum / (double)fN;
}

// Definition of the propagation over a time interval of n time steps,
// whose expectation value is computed in the Metropolis::ComputePropagatorEstimators.
// path: path on which the observable is evaluated
// n: time index taking values in the range of integers from 0 to fN -1.
double Metropolis::EvaluatePropagator(const std::vector<double>& path, int n) const
{
  double partial_sum = 0.0;
  if (path.size() != fN) {
    std::cout << "ERROR: wrong path size in EvaluatePropagator input\n";
    throw 1;
  }
  for (int k = 0; k < fN; k++) partial_sum += path[k] * path[(k + n) % fN];
  return partial_sum / (double)fN;
}

#endif
