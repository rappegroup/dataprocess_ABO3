#ifndef atom_h
#define atom_h
#include <iostream>
#include <list>
#include "polarconfig.h"
typedef struct Atom{
	double position[3];
	char type;
	double charge;
	int tick;
}atom;
void sort(double* input,int dim);
double* distance(atom *a,atom *b,double *p);
double norm(double* p,int dim);
/*no periodical boudary condition*/
double* distance(double*a ,double* b);
double far(atom* a,atom* b,double* p);
double* displace_average_A(atom *A,atom *oxygen,double* p,int cell);
double* displace_average_B(atom *B,atom *oxygen,double* p,int cell);
double* displace_average_Ba(atom *A,atom *oxygen,double* p,int cell);
double* displace_average_Ca(atom *A,atom *oxygen,double* p,int cell);
double* polar_average(atom *A,atom *B,atom *oxygen,double* p,int cell);
double displace_average_B_scalar(atom *B,atom *oxygen,double* p,int cell);
double displace_average_Ba_scalar(atom *A,atom *oxygen,double* p,int cell);
double displace_average_Ca_scalar(atom *A,atom *oxygen,double* p,int cell);
void sum_together(double*,double*,int);
int* changeindex(int index,int cell);
int  changeback(int x,int y,int z,int cell);
double average(std::list<double> &input);
double variance(std::list<double> &input);
int* neighbor_o_forA(int index,int cell);
int* neighbor_o_forB(int index,int cell);
int* neighbor_A_forB(int index,int cell);
double tiltangle(atom* a, atom* b,atom* c,double* p);
void analyzepolar(atom* A,atom* B,atom* oxygen, double* period, int cell);
void outpolar();
#endif
