////////////////////////////////////////////////////////////////////////
/// \file my4Vector.h
/// \brief Header file for the definition of the class my4Vector
///
/// Header file containing the definitions of the attributes and members
/// of the class my4Vector. Further comments may be found in the
/// implementation file of this class.
////////////////////////////////////////////////////////////////////////
#ifndef MY4VECTOR_H
#define MY4VECTOR_H

#include <vector>

/// my4Vector class
///
/// This class represents the position 4-vector on the 4D lattice.
/// It includes the methods my4Vector::Move and my4Vector::Offset, which
/// allow to move with respect to a given point with built-in periodic
/// boundary conditions.
class my4Vector
{
 private:
  std::vector<int>
      fPosition;  ///< Vector containing the position 4-vector in terms of integer indeces
  std::vector<int> fNCells;  ///< Vector containing the lattice dimensions in the 4 dimensions

  /// Set
  ///
  /// Set a single position in a given component mu to pos.
  /// \param pos new position as integer index
  /// \param mu selected component
  void Set(int pos, int mu);

 public:
  /// Constructor
  ///
  /// \param position vector containing the position 4-vector in terms of integer indeces
  /// \param ncells vector containing the lattice dimensions in the 4 dimensions
  my4Vector(std::vector<int> position, std::vector<int> ncells);
  my4Vector() = delete;

  /// Destructor
  ~my4Vector();

  /// [] overloading
  ///
  /// Overloading of the [] operator to access the value of position index.
  /// \param mu component whose value is accessed
  /// \return position index in the mu component
  int operator[](int mu) const;

  /// \return fPosition
  std::vector<int> GetPosition() const;

  /// \return fNCells
  std::vector<int> GetNCells() const;

  /// Move
  ///
  /// Move the current position of nsteps in the mu component. This function
  /// implements the cyclic boundary conditions on both time and space coordinates.
  /// \param nsteps number of integer steps to move from current position
  /// \param mu component to move
  void Move(int nsteps, int mu);

  /// Offset
  ///
  /// Create a copy of the current position with an offset of nsteps in the mu component. This
  /// method implements the cyclic boundary conditions on both time and space coordinates. \param
  /// nsteps number of integer steps to move from current position \param mu component to move
  my4Vector Offset(int nsteps, int mu) const;
};

#endif
