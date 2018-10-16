#include "atom.h"
#include <iostream>
#include <fstream>
#include <string>
#include <new>
#include <algorithm>
#include <map>
#include <cmath>
#include <vector>
#include "space.h"
int main(int argc,char** argv){
	//caculating the displacement for ABO_x3
	int cell=std::stoi(argv[1]);
	std::fstream dump;
	std::fstream result;
	std::fstream calist;
	dump.open(argv[2],std::fstream::in);
	calist.open(argv[3],std::fstream::in);
	atom* A=new atom[cell*cell*cell];
	double* dispba;
	double* dispca;
	double disp_scalar;
	atom* B=new atom[cell*cell*cell];
	double* dispB;
	double* polar;
	atom* oxygen=new atom[3*cell*cell*cell];
	for(size_t i=0;i<cell*cell*cell;i++){
		A[i].type='b';
		A[i].charge=2.70;
	}
	int ca_num;
	while(calist>>ca_num){
		A[ca_num].type='c';
	}
	for(size_t i=0;i<cell*cell*cell;i++){
		B[i].type='z';
		B[i].charge=6.30;
	}
	for(size_t i=0;i<3*cell*cell*cell;i++){
		oxygen[i].type='o';
		oxygen[i].charge=-3.00;
	}
	std::string la_pattern="ITEM: BOX BOUNDS pp pp pp";
	std::string coord_pattern="ITEM: ATOMS x y z ";
	std::list<double> disp_allba_x;
	std::list<double> disp_allba_y;
	std::list<double> disp_allba_z;
	std::list<double> disp_allca_x;
	std::list<double> disp_allca_y;
	std::list<double> disp_allca_z;
	std::list<double> disp_ca_scalar;
	std::list<double> disp_ba_scalar;
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
	int* index;
	int a;
	int b;
	int c;
	double angle;
	double period[3]={0,0,0};
	double x1,x2;
	size_t signal=0;
	for(std::string line;getline(dump,line);){
		if(la_pattern==line){
			for(size_t i=0;i<3;i++){
			dump>>x1;
			dump>>x2;
			period[i]=x2-x1;
			}
			la_x.push_back(period[0]/cell);
			la_y.push_back(period[1]/cell);
			la_z.push_back(period[2]/cell);
		}
	  if(coord_pattern==line || line.find(coord_pattern)!=std::string::npos){
			for(size_t i=0;i<cell*cell*cell;i++){
				for(size_t j=0;j<3;j++){
					dump>>A[i].position[j];
				}
			}
			for(size_t i=0;i<cell*cell*cell;i++){
				for(size_t j=0;j<3;j++){
					dump>>B[i].position[j];
				}
			}
			for(size_t i=0;i<3*cell*cell*cell;i++){
				for(size_t j=0;j<3;j++){
				dump>>oxygen[i].position[j];
				}
			}
            std::cout<<Pentahedron(A,oxygen,0,cell,period)<<std::endl;
			dispB=displace_average_B(B,oxygen,period,cell);
			disp_scalar=displace_average_B_scalar(B,oxygen,period,cell);
			disp_B_scalar.push_back(disp_scalar);
			disp_allB_x.push_back(dispB[0]);
			disp_allB_y.push_back(dispB[1]);
			disp_allB_z.push_back(dispB[2]);
			dispba=displace_average_Ba(A,oxygen,period,cell);
			disp_scalar=displace_average_Ba_scalar(A,oxygen,period,cell);
			disp_ba_scalar.push_back(disp_scalar);
			disp_allba_x.push_back(dispba[0]);
			disp_allba_y.push_back(dispba[1]);
			disp_allba_z.push_back(dispba[2]);
			dispca=displace_average_Ca(A,oxygen,period,cell);
			disp_scalar=displace_average_Ca_scalar(A,oxygen,period,cell);
			disp_ca_scalar.push_back(disp_scalar);
			disp_allca_x.push_back(dispca[0]);
			disp_allca_y.push_back(dispca[1]);
			disp_allca_z.push_back(dispca[2]);
			polar=polar_average(A,B,oxygen,period,cell);
			//sort(polar,3);
			px.push_back(polar[0]);
			py.push_back(polar[1]);
			pz.push_back(polar[2]);
			//compute the tilt angle now;
			for(size_t i=0;i<cell*cell*cell;i++){
				index=changeindex(i,cell);
				a=changeback(index[0],index[1],index[2],cell);
				b=changeback(index[0],index[1],index[2]+1,cell);
				c=changeback(index[0],index[1],index[2]+2,cell);
				angle=tiltangle(a+oxygen,b+oxygen,c+oxygen,period);
				tilt_angle.push_back(angle);
				tilt_angle_one.push_back(angle);
			}
			for(size_t i=0;i<cell*cell*cell;i++){
				index=changeindex(i,cell);
				a=changeback(index[0],index[1],index[2],cell);
				b=changeback(index[0],index[1]+1,index[2],cell);
				c=changeback(index[0],index[1]+2,index[2],cell);
				angle=tiltangle(cell*cell*cell+a+oxygen,cell*cell*cell+b+oxygen,c+oxygen+cell*cell*cell,period);
				tilt_angle.push_back(angle);
				tilt_angle_two.push_back(angle);
			}
		for(size_t i=0;i<cell*cell*cell;i++){
				index=changeindex(i,cell);
				a=changeback(index[0],index[1],index[2],cell);
				b=changeback(index[0]+1,index[1],index[2],cell);
				c=changeback(index[0]+2,index[1],index[2],cell);
				angle=tiltangle(2*cell*cell*cell+a+oxygen,2*cell*cell*cell+b+oxygen,c+oxygen+2*cell*cell*cell,period);
				tilt_angle.push_back(angle);
				tilt_angle_three.push_back(angle);
			}
		}
		}
	std::fstream fileout;
	fileout.open("result.txt",std::fstream::out);
	fileout<<"the average lattice constant is:"<<std::endl;
	fileout<<average(la_x)<<" "<<average(la_y)<<" "<<average(la_z)<<std::endl;
	fileout<<"the average B site cations displacement is:"<<std::endl;
	fileout<<fabs(average(disp_allB_x))<<" "<<fabs(average(disp_allB_y))<<" "<<fabs(average(disp_allB_z))<<std::endl;
	fileout<<"the average Asite one displacement is:"<<std::endl;
	fileout<<fabs(average(disp_allba_x))<<" "<<fabs(average(disp_allba_y))<<" "<<fabs(average(disp_allba_z))<<std::endl;
	fileout<<"the average Asite two displacement is:"<<std::endl;
	fileout<<fabs(average(disp_allca_x))<<" "<<fabs(average(disp_allca_y))<<" "<<fabs(average(disp_allca_z))<<std::endl;
	fileout<<"the average O6 tilt angle is:"<<std::endl;
	fileout<<(fabs(average(tilt_angle)))<<std::endl;
	fileout<<"the averaget O6 tilt angle in three direction is:"<<std::endl;
	fileout<<(fabs(average(tilt_angle_one)))<<" "<<fabs(average(tilt_angle_two))<<" "<<fabs(average(tilt_angle_three))<<std::endl;
	fileout<<"the scalar average B site displacement is:"<<std::endl;
	fileout<<average(disp_B_scalar)<<std::endl;
	fileout<<"the scalar average Asite one dispalcement is:"<<std::endl;
	fileout<<average(disp_ba_scalar)<<std::endl;
	fileout<<"the scalar average Asite two displacement is:"<<std::endl;
	fileout<<average(disp_ca_scalar)<<std::endl;
	std::vector<double> pall(3,0.0);
	std::vector<double> var(3,0.0);
	pall[0]=std::fabs(average(px));
	pall[1]=std::fabs(average(py));
	pall[2]=std::fabs(average(pz));
	std::fstream pout;
	pout.open("polar.txt",std::fstream::out);
	std::list<double>::iterator pyi=py.begin();
	std::list<double>::iterator pzi=pz.begin();
	for(std::list<double>::iterator pxi=px.begin();pxi!=px.end();pxi++){
		pout<<*(pxi)<<" "<<*pyi<<" "<<*pzi<<std::endl;
		pyi++;
		pzi++;
	}
	var[0]=variance(px);
	var[1]=variance(py);
	var[2]=variance(pz);
	std::map <double,double> good;
	for(size_t i=0;i<3;i++){
		good.insert(good.end(),std::pair <double,double> (pall[i],var[i]));
	}
//	sort(pall.begin(),pall.end());
	fileout<<"the polarization is (absolute value):"<<std::endl;
	fileout<<pall[0]<<" "<<pall[1]<<" "<<pall[2]<<std::endl;
	fileout<<"polarization variance is:"<<std::endl;
	fileout<<good[pall[0]]<<" "<<good[pall[1]]<<" "<<good[pall[2]]<<std::endl;
	dump.close();
	calist.close();
	fileout.close();
	return 0;
}
