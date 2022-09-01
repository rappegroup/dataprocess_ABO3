#ifndef polarconfig_h
#define polarconfig_h
#include <list>
#include <vector>
#include <iostream>
#include <queue>
namespace polarconfig
{
   extern std::list<double> disp_B_scalar;
   extern std::list<double> disp_allB_x;
   extern std::list<double> disp_allB_y;
   extern std::list<double> disp_allB_z;
   extern std::list<double> tilt_angle;
   extern std::list<double> tilt_angle_one;
   extern std::list<double> tilt_angle_two;
   extern std::list<double> tilt_angle_three;
   extern std::list<double> la_x;
   extern std::list<double> la_y;
   extern std::list<double> la_z;
   extern std::list<double> px;
   extern std::list<double> py;
   extern std::list<double> pz;
   extern int** mapunit;
   extern int* map1D;
   extern int** mapunitA;
   extern int* map1DA;
   extern std::vector<std::list<double> > px_local;
   extern std::vector<std::list<double> > py_local;
   extern std::vector<std::list<double> > pz_local;
   extern std::vector<double> epsilon_x;
   extern std::vector<double> epsilon_y;
   extern std::vector<double> epsilon_z;
   extern double tilt_one_ave;
   extern double tilt_two_ave;
   extern double tilt_three_ave;
   extern double temperature;
   extern int Nx;
   extern int Ny;
   extern int Nz;
   extern int cell;
   extern int steps;
}
#endif
