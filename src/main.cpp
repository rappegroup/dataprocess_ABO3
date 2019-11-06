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
#include <list>
#include <mpi.h>
#include <sstream>
#include <cstdlib>
int main(){
	//caculating the displacement for ABO_x3
  MPI_Init(NULL,NULL);
  int& cell=polarconfig::cell;
  std::fstream dump;
  std::fstream calist;
  std::string dumpfile;
  int velocity_on=1;
  int polarization_on=1;
  int position_variance_on=1;
  int local_die=0;
  std::string calistfile;
  int world_rank;
  int world_size;
  MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
  std::list<double*> ve_list;
  if(world_rank==0){
    info(cell,dumpfile,calistfile,velocity_on,polarization_on,polarconfig::temperature,position_variance_on,local_die);
	std::cout<<"the temperature now is: "<<polarconfig::temperature<<std::endl;
	}
  std::cout<<"I am here fine:"<<std::endl;
  MPI_Bcast(&cell,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&polarconfig::temperature,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&velocity_on,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&polarization_on,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&position_variance_on,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&local_die,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Comm_size(MPI_COMM_WORLD,&world_size);
  int MPI_LOOP_COUNT=ceil((cell*cell*cell+0.0)/world_size);
  double* ve_temp;
	size_t v_count=0;
  if(world_rank==0){
	dump.open(dumpfile.c_str(),std::fstream::in);
	calist.open(calistfile.c_str(),std::fstream::in);
  }
	atom* A=new atom[cell*cell*cell];
	atom* B=new atom[cell*cell*cell];
	atom* oxygen=new atom[3*cell*cell*cell];
  clock_t begin=clock();
	for(size_t i=0;i<cell*cell*cell;i++){
		A[i].type='b';
		A[i].charge[0]=2.90;
    A[i].charge[1]=2.90;
    A[i].charge[2]=2.90;
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
  atom atom_demo;
  int blockcounts[4]={3,3,1,1};
  MPI_Datatype types[4];
  MPI_Aint displs[4];
  MPI_Datatype MPI_atom;
  MPI_Address(atom_demo.position,&displs[0]);
  MPI_Address(atom_demo.charge,&displs[1]);
  MPI_Address(&atom_demo.type,&displs[2]);
  MPI_Address(&atom_demo.tick,&displs[3]);
  for(int i=3;i>=0;i--){
    displs[i]=displs[i]-displs[0];
  }
  types[0]=MPI_DOUBLE;
  types[1]=MPI_DOUBLE;
  types[2]=MPI_CHAR;
  types[3]=MPI_INT;
  MPI_Type_struct(4,blockcounts,displs,types,&MPI_atom);
  MPI_Type_commit(&MPI_atom);
	std::string la_pattern="ITEM: BOX BOUNDS pp pp pp";
	std::string coord_pattern="ITEM: ATOMS x y z ";
	double period[3]={0,0,0};
	double x1,x2;
	size_t signal=0;
  int read_success;
  int getnewframe=0;
  size_t read_bound=60000;
  std::string line="0";
	do
   {
   if(world_rank==0){
    if(getline(dump,line)){
      read_success=1;
    }
    else{
      read_success=0;
     }
     }
   else{
    }
   MPI_Bcast(&read_success,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
   if(signal==read_bound){
    read_success=0;
    MPI_Bcast(&read_success,1,MPI_INT,0,MPI_COMM_WORLD);
   }
   if(read_success==1){
    //continue working on it;
   }
   else{
    break;
   }
   if(world_rank==0){
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
      getnewframe=1;
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
      }
      }
      else{
      //doing nothing, waiting for root processor finish reading.
      };
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Bcast(&getnewframe,1,MPI_INT,0,MPI_COMM_WORLD);
      if(getnewframe==1){
      MPI_Bcast(period,3,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(A,cell*cell*cell,MPI_atom,0,MPI_COMM_WORLD);
      MPI_Bcast(B,cell*cell*cell,MPI_atom,0,MPI_COMM_WORLD);
      MPI_Bcast(oxygen,3*cell*cell*cell,MPI_atom,0,MPI_COMM_WORLD);
      MPI_Barrier(MPI_COMM_WORLD);
			if(polarization_on){
				analyzepolar(A,B,oxygen,period,cell);
			}
      if(position_variance_on){
        analyzeposition_variance(A,B,oxygen,period,cell,signal);
      }
      getnewframe=0;
      }
      else{
      }
	}while(true);
  clock_t end=clock();
  double use_secs = double(end - begin) / CLOCKS_PER_SEC;
  std::cout<<"The total time spend is: "<<use_secs<<std::endl;
  MPI_Barrier(MPI_COMM_WORLD);
	if(polarization_on){
        if(world_rank==0){
	     	outpolar();
    	}
	}
	if(velocity_on){
		autospeed(ve_list,cell);
	}
	dump.close();
	calist.close();
  MPI_Barrier(MPI_COMM_WORLD);
  if(local_die){
  double lx=average(polarconfig::la_x);
  double ly=average(polarconfig::la_y);
  double lz=average(polarconfig::la_z);
  double volume=lx*ly*lz;
  MPI_Bcast(&volume,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  calculate_local_die(cell,volume,polarconfig::temperature);
  }
  if(position_variance_on){
  calculate_local_variance(cell,polarconfig::temperature);
  }
  MPI_Finalize();
	return 0;
}
