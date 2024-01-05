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
#include "searchunit.h"
int main(){
	//caculating the displacement for ABO_x3
  MPI_Init(NULL,NULL);
  int& Nx=polarconfig::Nx;
  int& Ny=polarconfig::Ny;
  int& Nz=polarconfig::Nz;
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
  std::string datafile="STO.data";
  MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
  std::list<double*> ve_list;
  if(world_rank==0){
    info(Nx,Ny,Nz,dumpfile,calistfile,velocity_on,polarization_on,polarconfig::temperature,position_variance_on,local_die);
	std::cout<<"the temperature now is: "<<polarconfig::temperature<<std::endl;
	}
  MPI_Bcast(&Nx,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&Ny,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&Nz,1,MPI_INT,0,MPI_COMM_WORLD);
  searchneighbor(datafile,Nx,Ny,Nz);
  if(world_rank==1){
  for(size_t i=0;i<Nx*Ny*Nz;i++){
    for(size_t j=0;j<6+8;j++){
      std::cout<<polarconfig::mapunit[i][j]<<" ";
    }
    std::cout<<std::endl;
  }
  std::cout<<" finished by searching "<<std::endl;
  std::cout<<" the periodical length is: "<<std::endl;
  std::cout<<Nx<<" "<<Ny<<" "<<Nz<<std::endl;
  std::cout<<"the Asite neighbors are: "<<std::endl;
  for(size_t i=0;i<Nx*Ny*Nz;i++){
    for(size_t j=0;j<12;j++){
      std::cout<<polarconfig::mapunitA[i][j]<<" ";
    }
    std::cout<<std::endl;
  }
  std::cout<<" finished by searching "<<std::endl;
  std::cout<<" the periodical length is: "<<std::endl;
  std::cout<<Nx<<" "<<Ny<<" "<<Nz<<std::endl;
  }
  MPI_Bcast(&polarconfig::temperature,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&velocity_on,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&polarization_on,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&position_variance_on,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&local_die,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Comm_size(MPI_COMM_WORLD,&world_size);
  int MPI_LOOP_COUNT=ceil((Nx*Ny*Nz+0.0)/world_size);
  double* ve_temp;
	size_t v_count=0;
  if(world_rank==0){
	dump.open(dumpfile.c_str(),std::fstream::in);
	calist.open(calistfile.c_str(),std::fstream::in);
  }
	atom* A=new atom[Nx*Ny*Nz];
	atom* B=new atom[Nx*Ny*Nz];
	atom* oxygen=new atom[3*Nx*Ny*Nz];
if( world_rank == 0){
    std::fstream chargefile;
    chargefile.open("CHARGE.dat", std::fstream::in);
    for(size_t i = 0; i < Nx * Ny * Nz; i++){
        A[i].type = 'a';
        for(size_t j = 0; j < 3; j++){
            chargefile >> A[i].charge[j];
        }
    }
    for(size_t i = 0; i < Nx * Ny * Nz; i++){
        B[i].type = 'b';
        for(size_t j = 0; j < 3; j++){
            chargefile >> B[i].charge[j];
        }
    }
    for(size_t i = 0; i < 3 * Nx * Ny * Nz; i++){
        oxygen[i].type = 'o';
        for(size_t j = 0; j < 3; j++){
            chargefile >> oxygen[i].charge[j];
        }
    }
}
      MPI_Barrier(MPI_COMM_WORLD);
  clock_t begin=clock();
  atom atom_demo;
  int blockcounts[4]={3,3,1,1};
  MPI_Datatype types[4];
  MPI_Aint displs[4];
  MPI_Datatype MPI_atom;
  MPI_Get_address(atom_demo.position,&displs[0]);
  MPI_Get_address(atom_demo.charge,&displs[1]);
  MPI_Get_address(&atom_demo.type,&displs[2]);
  MPI_Get_address(&atom_demo.tick,&displs[3]);
  for(int i=3;i>=0;i--){
    displs[i]=displs[i]-displs[0];
  }
  types[0]=MPI_DOUBLE;
  types[1]=MPI_DOUBLE;
  types[2]=MPI_CHAR;
  types[3]=MPI_INT;
  MPI_Type_create_struct(4,blockcounts,displs,types,&MPI_atom);
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
			polarconfig::la_x.push_back(period[0]/Nx);
			polarconfig::la_y.push_back(period[1]/Ny);
			polarconfig::la_z.push_back(period[2]/Nz);
			if(velocity_on){
				ve_temp=new double [Nx*Ny*Nz*5*3];
				ve_list.push_back(ve_temp);
				v_count=0;
			}
		}
	  if(coord_pattern==line || line.find(coord_pattern)!=std::string::npos){
      getnewframe=1;
			for(size_t i=0;i<Nx*Ny*Nz;i++){
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
			for(size_t i=0;i<Nx*Ny*Nz;i++){
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
			for(size_t i=0;i<3*Nx*Ny*Nz;i++){
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
      MPI_Bcast(A,Nx*Ny*Nz,MPI_atom,0,MPI_COMM_WORLD);
      MPI_Bcast(B,Nx*Ny*Nz,MPI_atom,0,MPI_COMM_WORLD);
      MPI_Bcast(oxygen,3*Nx*Ny*Nz,MPI_atom,0,MPI_COMM_WORLD);
      MPI_Barrier(MPI_COMM_WORLD);
			if(polarization_on){
				polar_calculate_search(A,B,oxygen,period,Nx,Ny,Nz);
        dispA_calculate_search(A,B,oxygen,period,Nx,Ny,Nz);
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
  MPI_Finalize();
	return 0;
}
