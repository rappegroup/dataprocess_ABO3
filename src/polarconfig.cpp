#include "polarconfig.h"
#include <list>
#include <iostream>
#include <queue>
namespace polarconfig
{
 std::list<double> disp_allba_x;
 std::list<double> disp_allba_y;
 std::list<double> disp_allba_z;
 std::list<double> disp_allca_x;
 std::list<double> disp_allca_y;
 std::list<double> disp_allca_z;
 std::list<double> disp_ca_scalar;
 std::list<double> disp_ba_scalar;
 std::list<double> disp_bi_scalar;
 std::list<double> disp_ti_scalar;
 std::list<double> disp_mg_scalar;
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
 std::vector<std::list<double> > disp_Asite_x;
 std::vector<std::list<double> > disp_Asite_y;
 std::vector<std::list<double> > disp_Asite_z;
 std::vector<std::list<double> > disp_Bsite_x;
 std::vector<std::list<double> > disp_Bsite_y;
 std::vector<std::list<double> > disp_Bsite_z;
 std::vector<std::queue<double> > dipole_x;
 std::vector<std::queue<double> > dipole_y;
 std::vector<std::queue<double> > dipole_z;
 double temperature;
 int cell;
}
