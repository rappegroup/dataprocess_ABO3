#include "polarconfig.h"
#include <list>
#include <iostream>
#include <queue>
namespace polarconfig
{
 std::list<double> disp_B_scalar;
 std::list<double> disp_allB_x;
 std::list<double> disp_allB_y;
 std::list<double> disp_allB_z;
 std::list<double> tilt_angle;
 std::list<double> tilt_angle_one;
 std::list<double> tilt_angle_two;
 std::list<double> tilt_angle_three;
 std::list<double> la_x;
 std::list<double> la_y;
 std::list<double> la_z;
 std::list<double> px;
 std::list<double> py;
 std::list<double> pz;
 std::vector<std::list<double> > px_local;
 std::vector<std::list<double> > py_local;
 std::vector<std::list<double> > pz_local;
 std::vector<double> epsilon_x;
 std::vector<double> epsilon_y;
 std::vector<double> epsilon_z;
 double tilt_one_ave;
 double tilt_two_ave;
 double tilt_three_ave;
 double temperature;
 int cell;
}
