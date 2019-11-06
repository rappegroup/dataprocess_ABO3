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
  MPI_Barrier(MPI_COMM_WORLD);
  if(local_die){
  calculate_local_die(cell,4.02*4.02*4.02,polarconfig::temperature);
  }
//  if(position_variance_on){
//  calculate_local_variance(cell,polarconfig::temperature);
//  }
  MPI_Finalize();
	return 0;
}
