#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <string>
#include <list>
double average(std::list<double>& listA){
  double sum=0.0;
  int count=0;
  for(std::list<double>::iterator a=listA.begin();a!=listA.end();a++){
    sum=sum+*a;
    count=count+1;
  }
  return sum/count;
}
int main(){
  MPI_Init(NULL,NULL);
  int simulationtime=100000/200;
  int cell=20;
  int world_rank,world_size;
  MPI_Comm_size(MPI_COMM_WORLD,&world_size);
  MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
  MPI_File mpifile;
  MPI_Status status;
  MPI_Offset offset;
  std::string filename="local_polar.bin";
  double* localpolar=new double[3];
  std::list<double> px,py,pz;
  MPI_File_open(MPI_COMM_WORLD,filename.c_str(),MPI_MODE_RDONLY,MPI_INFO_NULL,&mpifile);
  double* polaraverage=new double [3*cell*cell*cell];
  for(size_t i=0;i<3*cell*cell*cell;i++){
    polaraverage[i]=0.0;
  }
  for(size_t loop=world_rank;loop<cell*cell*cell;loop=loop+world_size){
    px.clear();
    py.clear();
    pz.clear();
    for(size_t frame=0;frame<simulationtime;frame++){
      offset=(frame*(3*cell*cell*cell)+loop*3)*sizeof(double);
      MPI_File_read_at(mpifile,offset,localpolar,3,MPI::DOUBLE,&status);
      px.push_back(localpolar[0]);
      py.push_back(localpolar[1]);
      pz.push_back(localpolar[2]);
    }
    polaraverage[loop*3+0]=average(px);
    polaraverage[loop*3+1]=average(py);
    polaraverage[loop*3+2]=average(pz);
  }
  double* reducepolar=new double [3*cell*cell*cell];
  MPI_Allreduce(polaraverage,reducepolar,3*cell*cell*cell,MPI::DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_File_close(&mpifile);
  MPI_Finalize();
}
