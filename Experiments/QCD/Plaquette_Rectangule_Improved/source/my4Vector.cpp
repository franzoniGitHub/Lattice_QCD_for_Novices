#include <vector>
#include "my4Vector.h"

// Constructor
my4Vector::my4Vector(std::vector<int> position, std::vector<int> ncells)
    : fPosition(position), fNCells(ncells)
{
  if (position.size() != 4 || ncells.size() != 4) throw 1;
  bool cond1 = false, cond2 = false;
  for (int i = 0; i < 4; i++) cond1 = cond1 || (ncells[i] <= 0);
  for (int i = 0; i < 4; i++) cond2 = cond2 || ((position[i] < 0) || (position[i] >= ncells[i]));
  if (cond1 || cond2) throw 1;
}

// Destructor
my4Vector::~my4Vector()
{
}

// Overloading of [] operator. It allows to get position in the mu coordinate
int my4Vector::operator[](int mu) const
{
  if (mu < 0 || mu >= 4) throw 1;
  return fPosition[mu];
}

// Set position in the mu coordinate
void my4Vector::Set(int pos, int mu)
{
  if (mu < 0 || mu >= 4) throw 1;
  if (pos < 0 || pos >= fNCells[mu]) throw 1;
  fPosition[mu] = pos;
}

// Getters
std::vector<int> my4Vector::GetPosition() const
{
  return fPosition;
}
std::vector<int> my4Vector::GetNCells() const
{
  return fNCells;
}

// Move the position of n-steps in the selected direction. This function
// implements the cyclic boundary conditions on both time and space coordinates
void my4Vector::Move(int nsteps, int mu)
{
  if (mu < 0 || mu >= 4) throw 1;
  if ((fPosition[mu] + nsteps) >= 0) {
    Set((fPosition[mu] + nsteps) % fNCells[mu], mu);
  } else {
    int k = 1;
    while (k * fNCells[mu] + fPosition[mu] + nsteps < 0) k++;
    Set(k * fNCells[mu] + fPosition[mu] + nsteps, mu);
  }
}

// Offset: it returns a my4Vector whose position is nsteps apart in the mu
// direction. It is the same as move, but without affecting the current position
// of the calling instance.
my4Vector my4Vector::Offset(int nsteps, int mu) const
{
  if (mu < 0 || mu >= 4) throw 1;
  std::vector<int> newPosition = fPosition;
  if ((fPosition[mu] + nsteps) >= 0) {
    newPosition[mu] = (fPosition[mu] + nsteps) % fNCells[mu];
  } else {
    int k = 1;
    while (k * fNCells[mu] + fPosition[mu] + nsteps < 0) k++;
    newPosition[mu] = k * fNCells[mu] + fPosition[mu] + nsteps;
  }
  return my4Vector(newPosition, fNCells);
}
