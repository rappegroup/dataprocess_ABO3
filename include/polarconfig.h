#ifndef polarconfig_h
#define polarconfig_h
#include <list>
#include <iostream>
namespace polarconfig
{
   extern std::list<double> disp_allba_x;
   extern std::list<double> disp_allba_y;
   extern std::list<double> disp_allba_z;
   extern std::list<double> disp_allca_x;
   extern std::list<double> disp_allca_y;
   extern std::list<double> disp_allca_z;
   extern std::list<double> disp_ca_scalar;
   extern std::list<double> disp_ba_scalar;
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
	 extern double temperature;
	 extern int cell;
}
#endif
