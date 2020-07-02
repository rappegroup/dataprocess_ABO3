#ifndef searchunit_h
#define searchunit_h
#include <string>
#include "atom.h"
void searchneighbor(std::string file,int cell);
void searchneighbor(std::string file,int Nx,int Ny,int Nz);
void polar_calculate_search(atom* A,atom* B,atom* oxygen,double* p,int cell);
void polar_calculate_search(atom* A,atom* B,atom* oxygen,double* p,int Nx,int Ny,int Nz);
#endif
