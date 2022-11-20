#include <armadillo>
#include <complex>
#include <iostream>
#include <string>
#include <vector>
#include "my4Vector.h"
#include "Path.h"
using namespace arma;

// Constructor
Path::Path(std::vector<int> ncells) : fNCells(ncells)
{
  if (ncells.size() != 4) throw 1;
  cx_dmat Identity(3, 3, fill::eye);
  std::vector<cx_dmat> lorentz;
  std::vector<std::vector<cx_dmat>> z_space;
  std::vector<std::vector<std::vector<cx_dmat>>> y_space;
  std::vector<std::vector<std::vector<std::vector<cx_dmat>>>> x_space;
  for (int mu = 0; mu < 4; mu++) lorentz.push_back(Identity);
  for (int i3 = 0; i3 < ncells[3]; i3++) z_space.push_back(lorentz);
  for (int i2 = 0; i2 < ncells[2]; i2++) y_space.push_back(z_space);
  for (int i1 = 0; i1 < ncells[1]; i1++) x_space.push_back(y_space);
  for (int i0 = 0; i0 < ncells[0]; i0++) fPath.push_back(x_space);
}

// Default constructor
Path::Path() : fNCells({1, 1, 1, 1})
{
  cx_dmat Identity(3, 3, fill::eye);
  std::vector<cx_dmat> lorentz;
  std::vector<std::vector<cx_dmat>> z_space;
  std::vector<std::vector<std::vector<cx_dmat>>> y_space;
  std::vector<std::vector<std::vector<std::vector<cx_dmat>>>> x_space;
  lorentz.push_back(Identity);
  z_space.push_back(lorentz);
  y_space.push_back(z_space);
  x_space.push_back(y_space);
  fPath.push_back(x_space);
}

// Destructor
Path::~Path()
{
}

// Reshape the lattice to a new size
void Path::Reshape(std::vector<int> ncells)
{
  if (ncells.size() != 4) throw 1;
  cx_dmat Identity(3, 3, fill::eye);
  std::vector<std::vector<std::vector<std::vector<std::vector<cx_dmat>>>>> temporary;
  std::vector<cx_dmat> lorentz;
  std::vector<std::vector<cx_dmat>> z_space;
  std::vector<std::vector<std::vector<cx_dmat>>> y_space;
  std::vector<std::vector<std::vector<std::vector<cx_dmat>>>> x_space;
  for (int mu = 0; mu < 4; mu++) lorentz.push_back(Identity);
  for (int i3 = 0; i3 < ncells[3]; i3++) z_space.push_back(lorentz);
  for (int i2 = 0; i2 < ncells[2]; i2++) y_space.push_back(z_space);
  for (int i1 = 0; i1 < ncells[1]; i1++) x_space.push_back(y_space);
  for (int i0 = 0; i0 < ncells[0]; i0++) temporary.push_back(x_space);
  fPath = temporary;
  fNCells = ncells;
}

// Get lattice dimensions
std::vector<int> Path::GetNCells() const
{
  return fNCells;
}

// Access to a lattice element
cx_dmat Path::operator()(const my4Vector& position, int mu) const
{
  bool condition = false;
  for (int i = 0; i < 4; i++) condition = condition || (position.GetNCells()[i] != fNCells[i]);
  if (condition) throw 1;
  int n0 = position[0], n1 = position[1], n2 = position[2], n3 = position[3];
  return fPath[n0][n1][n2][n3][mu];
}

// Access/Set lattice element
cx_dmat& Path::Set(const my4Vector& position, int mu)
{
  bool condition = false;
  for (int i = 0; i < 4; i++) condition = condition || (position.GetNCells()[i] != fNCells[i]);
  if (condition) throw 1;
  int n0 = position[0], n1 = position[1], n2 = position[2], n3 = position[3];
  return fPath[n0][n1][n2][n3][mu];
}

// Print path on screen: it prints the matrix determinant, too, as a check
void Path::Print() const
{
  std::cout << "\n\n#########################################################\n";
  for (int i0 = 0; i0 < fNCells[0]; i0++) {
    for (int i1 = 0; i1 < fNCells[1]; i1++) {
      for (int i2 = 0; i2 < fNCells[2]; i2++) {
        for (int i3 = 0; i3 < fNCells[3]; i3++) {
          for (int mu = 0; mu < 4; mu++) {
            std::cout << "***************\n";
            std::cout << "{" << i0 << ", " << i1 << ", " << i2 << ", " << i3 << ", " << mu << "}\n";
            fPath[i0][i1][i2][i3][mu].print();
            std::cout << det(fPath[i0][i1][i2][i3][mu]) << std::endl;
            std::cout << endl;
          }
        }
      }
    }
  }
  std::cout << "#########################################################\n\n";
}