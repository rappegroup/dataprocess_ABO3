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
#include <math.h>
void Global_polar(double* changeP,atom* AI,atom* A0,atom* BI,atom* B0,atom* OxyI,atom* Oxy0,double* period,int cell){
    int world_rank;
    int world_size;
    double temp;
    double dp[3]={0.0,0.0,0.0};
    MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
    MPI_Comm_size(MPI_COMM_WORLD,&world_size);
    for(size_t i=world_rank;i<cell*cell*cell;i=i+world_size){
        for(size_t j=0;j<3;j++){
            temp=AI[i].position[j]-A0[i].position[j];
            temp=temp-round(temp/period[j])*period[j];
            dp[j]=dp[j]+(temp)*AI[i].charge[j];
            temp=BI[i].position[j]-B0[i].position[j];
            temp=temp-round(temp/period[j])*period[j];
            dp[j]=dp[j]+(temp)*BI[i].charge[j];
        }
    }
    for(size_t i=world_rank;i<3*cell*cell*cell;i=i+world_size){
        for(size_t j=0;j<3;j++){
            temp=OxyI[i].position[j]-Oxy0[i].position[j];
            temp=temp-round(temp/period[j])*period[j];
            dp[j]=dp[j]+(temp)*OxyI[i].charge[j];
        }
    }
    for(size_t i=0;i<3;i++){
        dp[i]=dp[i]/period[0]/period[1]/period[2]*16.02;
    }
    double dpreduce[3]={0.0,0.0,0.0};
    MPI_Allreduce(dp,dpreduce,3,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    for(size_t i=0;i<3;i++){
        changeP[i]=dpreduce[i];
    }
}
int main(){
	//caculating the displacement for ABO_x3
  MPI_Init(NULL,NULL);
  int& cell=polarconfig::cell;
  std::fstream calist;
  std::string dumpfile;
  std::fstream chargefile;
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
    info(cell,polarconfig::Nx,polarconfig::Ny,polarconfig::Nz,dumpfile,calistfile,velocity_on,polarization_on,polarconfig::temperature,position_variance_on,local_die,polarconfig::steps);
	std::cout<<"the temperature now is: "<<polarconfig::temperature<<" "<<polarconfig::Nx<<" "<<polarconfig::Ny<<" "<<polarconfig::Nz<<" "<<polarconfig::steps<<std::endl;
	}
  MPI_Bcast(&cell,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&polarconfig::Nz,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&polarconfig::Ny,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&polarconfig::Nx,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&polarconfig::temperature,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&velocity_on,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&polarization_on,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&position_variance_on,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&local_die,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&polarconfig::steps,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Comm_size(MPI_COMM_WORLD,&world_size);
  int MPI_LOOP_COUNT=ceil((cell*cell*cell+0.0)/world_size);
  double* ve_temp;
	size_t v_count=0;
  if(world_rank==0){
	calist.open(calistfile.c_str(),std::fstream::in);
  }
	atom** A=new atom* [polarconfig::steps];
    double** Ashadow=new double* [polarconfig::steps];
	atom** B=new atom* [polarconfig::steps];
    double** Bshadow=new double* [polarconfig::steps];
	atom** oxygen=new atom* [polarconfig::steps];
    double** polarshadow=new double* [polarconfig::steps];
    double** period=new double* [polarconfig::steps];
    double** globalpolar=new double* [polarconfig::steps];
    for(size_t i=0;i<polarconfig::steps;i++){
        A[i]=new atom [polarconfig::Nx*polarconfig::Ny*polarconfig::Nz];
        Ashadow[i]=new double [3*polarconfig::Nx*polarconfig::Ny*polarconfig::Nz];
        B[i]=new atom [polarconfig::Nx*polarconfig::Ny*polarconfig::Nz];
        Bshadow[i]=new double [3*polarconfig::Nx*polarconfig::Ny*polarconfig::Nz];
        oxygen[i]=new atom [3*polarconfig::Nx*polarconfig::Ny*polarconfig::Nz];
        polarshadow[i]=new double [3*polarconfig::Nx*polarconfig::Ny*polarconfig::Nz];
        period[i]=new double [3];
        globalpolar[i]=new double [3];
    }
    std::cout<<"I am here 2"<<std::endl;
  clock_t begin=clock();
  if(world_rank==0){
    chargefile.open("CHARGE.dat",std::fstream::in);
	for(size_t i=0;i<cell*cell*cell;i++){
		A[0][i].type='b';
        for(size_t j=0;j<3;j++){
        chargefile>>A[0][i].charge[j];
        }
  	}
	for(size_t i=0;i<cell*cell*cell;i++){
		B[0][i].type='z';
        for(size_t j=0;j<3;j++){
            chargefile>>B[0][i].charge[j];
            }
	}
	for(size_t i=0;i<3*cell*cell*cell;i++){
	oxygen[0][i].type='o';
    for(size_t j=0;j<3;j++){
        chargefile>>oxygen[0][i].charge[j];
    }
	}
    chargefile.close();
    for(size_t i=0;i<polarconfig::steps;i++){
        for(size_t j=0;j<cell*cell*cell;j++){
        for(size_t k=0;k<3;k++){
            A[i][j].charge[k]=A[0][j].charge[k];
            B[i][j].charge[k]=B[0][j].charge[k];
            }
        }
        std::cout<<"i="<<i<<std::endl;
    }
    for(size_t i=0;i<polarconfig::steps;i++){
        for(size_t j=0;j<3*cell*cell*cell;j++){
        for(size_t k=0;k<3;k++){
            oxygen[i][j].charge[k]=oxygen[0][j].charge[k];
            }
        }
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  std::cout<<"I am here 3 "<<world_rank<<std::endl;
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
 MPI_Barrier(MPI_COMM_WORLD);
 /*Start to Read MD*/
  std::cout<<"I am here"<< world_rank <<std::endl;
 if(world_rank==0){
 std::fstream dump;
 dump.open(dumpfile.c_str(),std::fstream::in);
 readMD(dump,polarconfig::Nx,polarconfig::Ny,polarconfig::Nz,period,A,B,oxygen,polarconfig::steps);
 }
 else{
 }
 MPI_Barrier(MPI_COMM_WORLD);
 std::cout<<"Finished Reading"<<std::endl;
 for(size_t i=0;i<polarconfig::steps;i++){
      MPI_Bcast(period[i],3,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(A[i],polarconfig::Nx*polarconfig::Ny*polarconfig::Nz,MPI_atom,0,MPI_COMM_WORLD);
      MPI_Bcast(B[i],polarconfig::Nx*polarconfig::Ny*polarconfig::Nz,MPI_atom,0,MPI_COMM_WORLD);
      MPI_Bcast(oxygen[i],3*polarconfig::Nx*polarconfig::Ny*polarconfig::Nz,MPI_atom,0,MPI_COMM_WORLD);
      if(world_rank==0){
      std::cout<<"Time Step:="<<i<<std::endl;
      }
      if(polarization_on){
      analyzepolar(A[i],B[i],oxygen[i],Ashadow[i],Bshadow[i],polarshadow[i],period[i],cell);
      displace_A_unit(A[i],oxygen[i],Ashadow[i],period[i],cell);
      displace_B_unit(B[i],oxygen[i],Bshadow[i],period[i],cell);
       }
      Global_polar(globalpolar[i],A[i],A[0],B[i],B[0],oxygen[i],oxygen[0],period[i],cell);
      if(polarization_on){
        if(world_rank==0){
            outpolar();
         }
     }
 }
 MPI_Barrier(MPI_COMM_WORLD);
 clock_t end=clock();
 double use_secs = double(end - begin) / CLOCKS_PER_SEC;
 std::cout<<"The total time spend is: "<<use_secs<<std::endl;
 MPI_File fpolar,fdispA,fdispB,Global_P;
 MPI_File_open(MPI_COMM_WORLD,"local_polar.bin",MPI_MODE_CREATE | MPI_MODE_WRONLY,MPI_INFO_NULL,&fpolar);
 MPI_File_open(MPI_COMM_WORLD,"polar_direction_A.bin",MPI_MODE_CREATE | MPI_MODE_WRONLY,MPI_INFO_NULL,&fdispA);
 MPI_File_open(MPI_COMM_WORLD,"polar_direction_B.bin",MPI_MODE_CREATE | MPI_MODE_WRONLY,MPI_INFO_NULL,&fdispB);
 MPI_File_open(MPI_COMM_WORLD,"Global_polar.bin",MPI_MODE_CREATE | MPI_MODE_WRONLY,MPI_INFO_NULL,&Global_P);
  MPI_Offset offset;
 MPI_Status status;
 for(size_t i=world_rank;i<polarconfig::steps;i=i+world_size){
  offset=i*(3*polarconfig::Nx*polarconfig::Ny*polarconfig::Nz)*sizeof(double);
  MPI_File_write_at_all(fpolar,offset,polarshadow[i],3*polarconfig::Nx*polarconfig::Ny*polarconfig::Nz,MPI_DOUBLE,&status);
  MPI_File_write_at_all(fdispA,offset,Ashadow[i],3*polarconfig::Nx*polarconfig::Ny*polarconfig::Nz,MPI_DOUBLE,&status);
  MPI_File_write_at_all(fdispB,offset,Bshadow[i],3*polarconfig::Nx*polarconfig::Ny*polarconfig::Nz,MPI_DOUBLE,&status);
  offset=i*3*sizeof(double);
  MPI_File_write_at_all(Global_P,offset,globalpolar[i],3,MPI_DOUBLE,&status);
 }
 MPI_File_close(&fpolar);
 MPI_File_close(&fdispA);
 MPI_File_close(&fdispB);
 MPI_File_close(&Global_P);
 clock_t end2=clock();
 use_secs = double(end2 - end) / CLOCKS_PER_SEC;
 std::cout<<"The IO time spend is: "<<use_secs<<std::endl;
 MPI_Barrier(MPI_COMM_WORLD);
 MPI_Finalize();
}
