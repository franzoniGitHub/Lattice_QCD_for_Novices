////////////////////////////////////////////////////////////////////////
/// \file Path.h
/// \brief Header file for the definition of the class Path
///
/// Header file containing the definitions of the attributes and members
/// of the class Path. Further comments may be found in the
/// implementation file of this class.
////////////////////////////////////////////////////////////////////////
#ifndef PATH_H
#define PATH_H

#include <armadillo>
#include <complex>
#include <vector>
#include "my4Vector.h"
using namespace arma;

/// Path class
///
/// This class represents a generic lattice configuration. The lattice considered
/// is a 4D space-time lattice, where at each position 4 3x3 complex matrices are located.
/// Periodic boundary conditions are built in thanks to the interface with my4Vector class.
/// Access and set methods are implemented so as to simplify the notations in
/// calculations of functions of the lattice sites.
class Path
{
 private:
  std::vector<std::vector<std::vector<std::vector<std::vector<cx_dmat>>>>>
      fPath;  ///< Lattice configuration stored as a 5D vector of Armadillo 3x3 complex matrices
  std::vector<int> fNCells;  ///< Vector containing the lattice dimensions in the 4 dimensions

 public:
  /// Path constructor
  ///
  /// Initialize with 3x3 identity matrices a lattice of dimensions specified in ncells.
  /// \param ncells vector containing the lattice dimensions in the 4 dimensions
  Path(std::vector<int> ncells);

  /// Default Path constructor
  ///
  /// Initialize with a 3x3 identity matrix a lattice with a single site.
  Path();

  /// Destructor
  ~Path();

  /// Reshape
  ///
  /// Reshape the Path instance with new lattice dimensions and assign 3x3 identity matrices in
  /// every site. Useful when the Path instance is built with the default constructor \param ncells
  /// vector containing the lattice dimensions in the 4 dimensions
  void Reshape(std::vector<int> ncells);

  /// \return fNCells
  std::vector<int> GetNCells() const;

  /// () overloading
  ///
  /// Overloading of the () operator to access to the complex matrix of a given position and
  /// polarization in the lattice \param x position 4-vector \param mu polarization direction
  /// \return 3x3 matrix in the given site and polarization of the lattice
  cx_dmat operator()(const my4Vector& x, int mu) const;

  /// Set
  ///
  /// Method to set/acces the complex matrix at a given position and polarization on the lattice. In
  /// order to use this method to set, its application should occur on the left of an assignment
  /// operator, using the returned reference: this aims at improving the code appearence in complex
  /// calculations. \param x position 4-vector \param mu polarization direction \return 3x3 matrix
  /// in the given site and polarization of the lattice
  cx_dmat& Set(const my4Vector& x, int mu);

  /// Print
  ///
  /// Print the lattice configuration on standard output. The matrix determinants are printed, too.
  /// Useful for debugging purposes.
  void Print() const;
};

#endif