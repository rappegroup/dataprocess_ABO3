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
#include "autospeed.h"
#include <queue>
int main(int argc,char** argv){
	//caculating the displacement for ABO_x3
	int& cell=polarconfig::cell;
	std::fstream dump;
	std::fstream result;
	std::fstream calist;
	std::string dumpfile;
	bool velocity_on=true;
	bool polarization_on=true;
	std::string calistfile;
	info(cell,dumpfile,calistfile,velocity_on,polarization_on,polarconfig::temperature);
	std::list<double*> ve_list;
	std::cout<<"the temperature now is: "<<polarconfig::temperature<<std::endl;
	double* ve_temp;
	size_t v_count=0;
	dump.open(dumpfile.c_str(),std::fstream::in);
	calist.open(calistfile.c_str(),std::fstream::in);
	atom* A=new atom[cell*cell*cell];
	atom* B=new atom[cell*cell*cell];
	atom* oxygen=new atom[3*cell*cell*cell];
	for(size_t i=0;i<cell*cell*cell;i++){
		A[i].type='b';
		A[i].charge[0]=2.90;
    A[i].charge[1]=2.90;
    A[i].charge[2]=2.90;
    polarconfig::disp_Asite_x.push_back(std::list<double>());
    polarconfig::disp_Asite_y.push_back(std::list<double>());
    polarconfig::disp_Asite_z.push_back(std::list<double>());
    polarconfig::disp_Bsite_x.push_back(std::list<double>());
    polarconfig::disp_Bsite_y.push_back(std::list<double>());
    polarconfig::disp_Bsite_z.push_back(std::list<double>());
    polarconfig::dipole_x.push_back(std::queue<double>());
    polarconfig::dipole_y.push_back(std::queue<double>());
    polarconfig::dipole_z.push_back(std::queue<double>());
	}
	int ca_num;
	while(calist>>ca_num){
		A[ca_num].type='c';
	}
	for(size_t i=0;i<cell*cell*cell;i++){
		B[i].type='z';
    for(size_t j=0;j<3;j++){
		B[i].charge[j]=6.70;
    }
	}
	for(size_t i=0;i<cell*cell*cell;i++){
		oxygen[i].type='o';
		oxygen[i].charge[0]=-2.40;
    oxygen[i].charge[1]=-2.40;
    oxygen[i].charge[2]=-4.80;
	}
  for(size_t i=cell*cell*cell;i<2*cell*cell*cell;i++){
    oxygen[i].type='o';
    oxygen[i].charge[0]=-2.40;
    oxygen[i].charge[1]=-4.80;
    oxygen[i].charge[2]=-2.40;
  }
  for(size_t i=2*cell*cell*cell;i<3*cell*cell*cell;i++){
    oxygen[i].type='o';
    oxygen[i].charge[0]=-4.80;
    oxygen[i].charge[1]=-2.40;
    oxygen[i].charge[2]=-2.40;
  }
	std::string la_pattern="ITEM: BOX BOUNDS pp pp pp";
	std::string coord_pattern="ITEM: ATOMS x y z ";
	double period[3]={0,0,0};
	double x1,x2;
	size_t signal=0;
	for(std::string line;getline(dump,line);){
		if(line.find(la_pattern)!=std::string::npos){
			std::cout<<signal++<<std::endl;
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
      //std::cout<<Pentahedron(A,oxygen,0,cell,period)<<std::endl;
			if(polarization_on){
				analyzepolar(A,B,oxygen,period,cell);
        displace_Asite_record(A,oxygen,period,cell);
        displace_Bsite_record(B,oxygen,period,cell);
			}
		}
	}
	if(polarization_on){
		outpolar();
	}
	if(velocity_on){
		autospeed(ve_list,cell);
	}
	dump.close();
	calist.close();
	return 0;
}
