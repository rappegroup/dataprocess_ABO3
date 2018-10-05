#ifndef space_h
#define space_h
#include "atom.h"
/*calculate the volume of a tetrahedron under periodical boudary condition*/
double tetrahedron(double** nei,int tickA,int tickO1,int tickO2,int tickO3,double* p);
double Pentahedron(atom* A,atom* oxygen,int tickA,int cell,double* p);
double* outprod(double* a,double* b);
double inner(double* a,double* b);
int convex(double** nei,double* p,int a,int b,int c);
int diagconvex(double** nei,double* p,int a,int b,int c);
#endif
