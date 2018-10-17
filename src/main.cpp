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
#include "interface.h"
#include "polarconfig.h"
int main(int argc,char** argv){
	//caculating the displacement for ABO_x3
	int cell;
	std::fstream dump;
	std::fstream result;
	std::fstream calist;
	std::string dumpfile;
	bool velocity_on=true;
	bool polarization_on=true;
	std::string calistfile;
	info(cell,dumpfile,calistfile,velocity_on,polarization_on);
	std::list<double*> ve_list;
	double* ve_temp;
	size_t v_count=0;
	dump.open(dumpfile,std::fstream::in);
	calist.open(calistfile,std::fstream::in);
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
			polarconfig::la_x.push_back(period[0]/cell);
			polarconfig::la_y.push_back(period[1]/cell);
			polarconfig::la_z.push_back(period[2]/cell);
			if(velocity_on){
				ve_temp=new double [cell*cell*cell*5*3];
				ve_list.push_back(ve_temp);
				v_count=0;
			}
		}
	  if(coord_pattern==line || line.find(coord_pattern)!=std::string::npos){
			for(size_t i=0;i<cell*cell*cell;i++){
				for(size_t j=0;j<3;j++){
					dump>>A[i].position[j];
					}
				for(size_t j=0;j<3;j++){
					if(velocity_on){
						dump>>ve_temp[v_count];
						v_count++;
					}
				}
				}
			for(size_t i=0;i<cell*cell*cell;i++){
				for(size_t j=0;j<3;j++){
					dump>>B[i].position[j];
				}
				for(size_t j=0;j<3;j++){
					if(velocity_on){
					dump>>ve_temp[v_count];
					v_count++;
					}
				}
			}
			for(size_t i=0;i<3*cell*cell*cell;i++){
				for(size_t j=0;j<3;j++){
				dump>>oxygen[i].position[j];
				}
				for(size_t j=0;j<3;j++){
				if(velocity_on){
					dump>>ve_temp[v_count];
					v_count++;
					}
				}
			}
		/*
		for(size_t i=0;i<cell*cell*cell;i++){
			std::cout<<A[i].position[0]<<" "<<A[i].position[1]<<" "<<A[i].position[2]<<std::endl;
		}
		for(size_t i=0;i<cell*cell*cell;i++){
			std::cout<<B[i].position[0]<<" "<<B[i].position[1]<<" "<<B[i].position[2]<<std::endl;
		}
		for(size_t i=0;i<3*cell*cell*cell;i++){
			std::cout<<oxygen[i].position[0]<<" "<<oxygen[i].position[1]<<" "<<oxygen[i].position[2]<<std::endl;
		}
		*/
      std::cout<<Pentahedron(A,oxygen,0,cell,period)<<std::endl;
			dispB=displace_average_B(B,oxygen,period,cell);
			disp_scalar=displace_average_B_scalar(B,oxygen,period,cell);
			polarconfig::disp_B_scalar.push_back(disp_scalar);
			polarconfig::disp_allB_x.push_back(dispB[0]);
			polarconfig::disp_allB_y.push_back(dispB[1]);
			polarconfig::disp_allB_z.push_back(dispB[2]);
			dispba=displace_average_Ba(A,oxygen,period,cell);
			disp_scalar=displace_average_Ba_scalar(A,oxygen,period,cell);
			polarconfig::disp_ba_scalar.push_back(disp_scalar);
			polarconfig::disp_allba_x.push_back(dispba[0]);
			polarconfig::disp_allba_y.push_back(dispba[1]);
			polarconfig::disp_allba_z.push_back(dispba[2]);
			dispca=displace_average_Ca(A,oxygen,period,cell);
			disp_scalar=displace_average_Ca_scalar(A,oxygen,period,cell);
			polarconfig::disp_ca_scalar.push_back(disp_scalar);
			polarconfig::disp_allca_x.push_back(dispca[0]);
			polarconfig::disp_allca_y.push_back(dispca[1]);
			polarconfig::disp_allca_z.push_back(dispca[2]);
			polar=polar_average(A,B,oxygen,period,cell);
			//sort(polar,3);
			polarconfig::px.push_back(polar[0]);
			polarconfig::py.push_back(polar[1]);
			polarconfig::pz.push_back(polar[2]);
			//compute the tilt angle now;
			for(size_t i=0;i<cell*cell*cell;i++){
				index=changeindex(i,cell);
				a=changeback(index[0],index[1],index[2],cell);
				b=changeback(index[0],index[1],index[2]+1,cell);
				c=changeback(index[0],index[1],index[2]+2,cell);
				angle=tiltangle(a+oxygen,b+oxygen,c+oxygen,period);
				polarconfig::tilt_angle.push_back(angle);
				polarconfig::tilt_angle_one.push_back(angle);
			}
			for(size_t i=0;i<cell*cell*cell;i++){
				index=changeindex(i,cell);
				a=changeback(index[0],index[1],index[2],cell);
				b=changeback(index[0],index[1]+1,index[2],cell);
				c=changeback(index[0],index[1]+2,index[2],cell);
				angle=tiltangle(cell*cell*cell+a+oxygen,cell*cell*cell+b+oxygen,c+oxygen+cell*cell*cell,period);
				polarconfig::tilt_angle.push_back(angle);
				polarconfig::tilt_angle_two.push_back(angle);
			}
		for(size_t i=0;i<cell*cell*cell;i++){
				index=changeindex(i,cell);
				a=changeback(index[0],index[1],index[2],cell);
				b=changeback(index[0]+1,index[1],index[2],cell);
				c=changeback(index[0]+2,index[1],index[2],cell);
				angle=tiltangle(2*cell*cell*cell+a+oxygen,2*cell*cell*cell+b+oxygen,c+oxygen+2*cell*cell*cell,period);
				polarconfig::tilt_angle.push_back(angle);
				polarconfig::tilt_angle_three.push_back(angle);
			}
		}
	}
	if(polarization_on){
		outpolar();
	}
	dump.close();
	calist.close();
	return 0;
}
