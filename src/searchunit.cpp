#include "atom.h"
#include "polarconfig.h"
#include <algorithm>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <mpi.h>
void searchneighbor(std::string file,int cell){
  polarconfig::mapunit=new int* [cell*cell*cell];
  polarconfig::map1D=new int [cell*cell*cell*(6+8)];
  for(size_t i=0;i<cell*cell*cell;i++){
    polarconfig::mapunit[i]=&polarconfig::map1D[i*(6+8)];//1 Fe with 6 Oxygen, 8 Bi;
  }
  MPI_Barrier(MPI_COMM_WORLD);
  int world_rank,world_size;
  MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
  MPI_Comm_size(MPI_COMM_WORLD,&world_size);
  if(world_rank==0){
  std::fstream fs;
  fs.open(file.c_str(),std::fstream::in);
  std::string temp;
  int tempint;
  double tempdouble;
  double tempdouble2;
  double period[3]={0.0,0.0,0.0};
  std::stringstream ss;
  double** asite=new double* [cell*cell*cell];
  double** bsite=new double* [cell*cell*cell];
  double** osite=new double* [cell*cell*cell*3];
  while(getline(fs,temp)){
  if(temp.find("xlo")!=std::string::npos){
    ss.str(temp);
    ss>>tempdouble;
    ss>>tempdouble2;
    period[0]=tempdouble2-tempdouble;
  }
  if(temp.find("ylo")!=std::string::npos){
    ss.str(temp);
    ss>>tempdouble;
    ss>>tempdouble2;
    period[1]=tempdouble2-tempdouble;
  }
  if(temp.find("zlo")!=std::string::npos){
    ss.str(temp);
    ss>>tempdouble;
    ss>>tempdouble2;
    period[2]=tempdouble2-tempdouble;
  }
  if(temp.find("Atoms")!=std::string::npos){
    getline(fs,temp);
    for(size_t i=0;i<cell*cell*cell;i++){
      getline(fs,temp);
      asite[i]=new double[3];
      ss.str(temp);
      for(size_t j=0;j<3;j++){
        ss>>tempint;
      }
      ss>>tempdouble;
      for(size_t j=0;j<3;j++){
        ss>>asite[i][j];
      }
      ss.clear();
    }
    for(size_t i=0;i<cell*cell*cell;i++){
      getline(fs,temp);
      bsite[i]=new double[3];
      ss.str(temp);
      for(size_t j=0;j<3;j++){
        ss>>tempint;
      }
      ss>>tempdouble;
      for(size_t j=0;j<3;j++){
        ss>>bsite[i][j];
      }
      ss.clear();
    }
    for(size_t i=0;i<3*cell*cell*cell;i++){
      getline(fs,temp);
      osite[i]=new double[3];
      ss.str(temp);
      for(size_t j=0;j<3;j++){
        ss>>tempint;
      }
      ss>>tempdouble;
      for(size_t j=0;j<3;j++){
        ss>>osite[i][j];
      }
      ss.clear();
    }
  }
  }
  /*search the neighbors using the sort*/
  std::vector<std::pair<double,int> > sequence(cell*cell*cell,std::pair<double,int>(0.0,0));
  for(size_t i=0;i<cell*cell*cell;i++){
    for(size_t j=0;j<cell*cell*cell;j++){
      tempdouble=far(bsite[i],asite[j],period);
      sequence[j]=std::pair<double,int>(tempdouble,j);
    }
    std::sort(sequence.begin(),sequence.end());
    for(size_t j=0;j<8;j++){
    polarconfig::mapunit[i][j]=sequence[j].second;
    }
  }
  std::vector<std::pair<double,int> > sequence_two(3*cell*cell*cell,std::pair<double,int>(0.0,0));
  for(size_t i=0;i<cell*cell*cell;i++){
    for(size_t j=0;j<cell*cell*cell*3;j++){
      tempdouble=far(bsite[i],osite[j],period);
      sequence_two[j]=std::pair<double,int>(tempdouble,j);
    }
    std::sort(sequence_two.begin(),sequence_two.end());
    for(size_t j=0;j<6;j++){
      polarconfig::mapunit[i][j+8]=sequence_two[j].second;
    }
  }
  }
  else{
  };
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(polarconfig::map1D,cell*cell*cell*(6+8),MPI_INT,0,MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
//  for(size_t i=0;i<cell*cell*cell;i++){
//    for(size_t j=8;j<8+6;j++){
//      std::cout<<polarconfig::mapunit[i][j]<<" ";
//    }
//    std::cout<<std::endl;
//  }
}
void searchneighbor(std::string file,int Nx,int Ny,int Nz){
  polarconfig::mapunit=new int* [Nx*Ny*Nz];
  polarconfig::map1D=new int [Nx*Ny*Nz*(6+8)];
  polarconfig::mapunitA=new int* [Nx*Ny*Nz];
  polarconfig::map1DA=new int [Nx*Ny*Nz*12];
  for(size_t i=0;i<Nx*Ny*Nz;i++){
    polarconfig::mapunit[i]=&polarconfig::map1D[i*(6+8)];//1 Fe with 6 Oxygen, 8 Bi;
    polarconfig::mapunitA[i]=&polarconfig::map1DA[i*12];
  }
  MPI_Barrier(MPI_COMM_WORLD);
  int world_rank,world_size;
  MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
  MPI_Comm_size(MPI_COMM_WORLD,&world_size);
  if(world_rank==0){
  std::fstream fs;
  fs.open(file.c_str(),std::fstream::in);
  std::string temp;
  int tempint;
  double tempdouble;
  double tempdouble2;
  double period[3]={0.0,0.0,0.0};
  std::stringstream ss;
  double** asite=new double* [Nx*Ny*Nz];
  double** bsite=new double* [Nx*Ny*Nz];
  double** osite=new double* [Nx*Ny*Nz*3];
  while(getline(fs,temp)){
  if(temp.find("xlo")!=std::string::npos){
    ss.str(temp);
    ss>>tempdouble;
    ss>>tempdouble2;
    period[0]=tempdouble2-tempdouble;
  }
  if(temp.find("ylo")!=std::string::npos){
    ss.str(temp);
    ss>>tempdouble;
    ss>>tempdouble2;
    period[1]=tempdouble2-tempdouble;
  }
  if(temp.find("zlo")!=std::string::npos){
    ss.str(temp);
    ss>>tempdouble;
    ss>>tempdouble2;
    period[2]=tempdouble2-tempdouble;
  }
  if(temp.find("Atoms")!=std::string::npos){
    getline(fs,temp);
    for(size_t i=0;i<Nx*Ny*Nz;i++){
      getline(fs,temp);
      asite[i]=new double[3];
      ss.str(temp);
      for(size_t j=0;j<3;j++){
        ss>>tempint;
      }
      ss>>tempdouble;
      for(size_t j=0;j<3;j++){
        ss>>asite[i][j];
      }
      ss.clear();
    }
    for(size_t i=0;i<Nx*Ny*Nz;i++){
      getline(fs,temp);
      bsite[i]=new double[3];
      ss.str(temp);
      for(size_t j=0;j<3;j++){
        ss>>tempint;
      }
      ss>>tempdouble;
      for(size_t j=0;j<3;j++){
        ss>>bsite[i][j];
      }
      ss.clear();
    }
    for(size_t i=0;i<3*Nx*Ny*Nz;i++){
      getline(fs,temp);
      osite[i]=new double[3];
      ss.str(temp);
      for(size_t j=0;j<3;j++){
        ss>>tempint;
      }
      ss>>tempdouble;
      for(size_t j=0;j<3;j++){
        ss>>osite[i][j];
      }
      ss.clear();
    }
  }
  }
  /*search the neighbors using the sort*/
  std::vector<std::pair<double,int> > sequence(Nx*Ny*Nz,std::pair<double,int>(0.0,0));
  for(size_t i=0;i<Nx*Ny*Nz;i++){
    for(size_t j=0;j<Nx*Ny*Nz;j++){
      tempdouble=far(bsite[i],asite[j],period);
      sequence[j]=std::pair<double,int>(tempdouble,j);
    }
    std::sort(sequence.begin(),sequence.end());
    for(size_t j=0;j<8;j++){
    polarconfig::mapunit[i][j]=sequence[j].second;
    }
  }
  std::vector<std::pair<double,int> > sequence_two(3*Nx*Ny*Nz,std::pair<double,int>(0.0,0));
  for(size_t i=0;i<Nx*Ny*Nz;i++){
    for(size_t j=0;j<Nx*Ny*Nz*3;j++){
      tempdouble=far(bsite[i],osite[j],period);
      sequence_two[j]=std::pair<double,int>(tempdouble,j);
    }
    std::sort(sequence_two.begin(),sequence_two.end());
    for(size_t j=0;j<6;j++){
      polarconfig::mapunit[i][j+8]=sequence_two[j].second;
    }
  }
  std::vector<std::pair<double,int> > sequence_three(3*Nx*Ny*Nz,std::pair<double,int>(0.0,0));
  for(size_t i=0;i<Nx*Ny*Nz;i++){
    for(size_t j=0;j<Nx*Ny*Nz*3;j++){
      tempdouble=far(asite[i],osite[j],period);
      sequence_three[j]=std::pair<double,int>(tempdouble,j);
    }
    std::sort(sequence_three.begin(),sequence_three.end());
    for(size_t j=0;j<12;j++){
    polarconfig::mapunitA[i][j]=sequence_three[j].second;
    }
  }
  }
  else{
  };
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(polarconfig::map1D,Nx*Ny*Nz*(6+8),MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(polarconfig::map1DA,Nx*Ny*Nz*12,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
}
void polar_calculate_search(atom *A,atom *B,atom *oxygen,double *p,int cell){
	double volume=1.0;
	for(size_t k=0;k<3;k++){
		volume=volume*(p[k]/cell);
	}
	double* dist;
	double* sum=new double[3];
  int world_rank;
  int world_size;
  MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
  MPI_Comm_size(MPI_COMM_WORLD,&world_size);
  int MPI_LOOP_COUNT=ceil((cell*cell*cell+0.0)/world_size);
  int i=0;
  MPI_File fh;
  MPI_File_open(MPI_COMM_WORLD,"local_polar.bin", MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_APPEND,MPI_INFO_NULL,&fh);
  MPI_Offset initial_offset;
  MPI_Offset offset;
  MPI_File_get_position(fh,&initial_offset);
  MPI_Status status;
  int* neighbor;
	for(size_t layer=0;layer<MPI_LOOP_COUNT;layer++){
    i=layer*world_size+world_rank;
    if(i<cell*cell*cell){
		for(size_t k=0;k<3;k++){
			sum[k]=0.0;
		}
		for(size_t j=0;j<6;j++){
			dist=distance(B+i,polarconfig::mapunit[i][j+8]+oxygen,p);
			for(size_t k=0;k<3;k++){
				dist[k]=dist[k]*((oxygen+polarconfig::mapunit[i][j+8])->charge[k])/2.0;

			}
			sum_together(sum,dist,3);
			delete [] dist;
		}
		for(size_t j=0;j<8;j++){
			dist=distance(B+i,A+polarconfig::mapunit[i][j],p);
			for(size_t k=0;k<3;k++){
				//now this guy turn into polar.
				dist[k]=dist[k]*((A+polarconfig::mapunit[i][j])->charge[k])/8.0;
			}
			sum_together(sum,dist,3);
			delete [] dist;
      }
    for(size_t m=0;m<3;m++){
      sum[m]=sum[m]/volume*16;
    }
    offset=3*sizeof(double)*i;
    MPI_File_write_at_all(fh,initial_offset+offset,sum,3,MPI_DOUBLE,&status);
    }
    MPI_Barrier(MPI_COMM_WORLD);
	}
  MPI_File_close(&fh);
}
void polar_calculate_search(atom *A,atom *B,atom *oxygen,double *p,int Nx,int Ny,int Nz){
	double volume=1.0;
	for(size_t k=0;k<3;k++){
		volume=volume*(p[k]);
	}
  volume=volume/Nx/Ny/Nz;
	double* dist;
	double* sum=new double[3];
  int world_rank;
  int world_size;
  MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
  MPI_Comm_size(MPI_COMM_WORLD,&world_size);
  int MPI_LOOP_COUNT=ceil((Nx*Ny*Nz+0.0)/world_size);
  int i=0;
  MPI_File fh;
  MPI_File_open(MPI_COMM_WORLD,"local_polar.bin", MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_APPEND,MPI_INFO_NULL,&fh);
  MPI_Offset initial_offset;
  MPI_Offset offset;
  MPI_File_get_position(fh,&initial_offset);
  MPI_Status status;
  int* neighbor;
	for(size_t layer=0;layer<MPI_LOOP_COUNT;layer++){
    i=layer*world_size+world_rank;
    if(i<Nx*Ny*Nz){
		for(size_t k=0;k<3;k++){
			sum[k]=0.0;
		}
		for(size_t j=0;j<6;j++){
			dist=distance(B+i,polarconfig::mapunit[i][j+8]+oxygen,p);
			for(size_t k=0;k<3;k++){
				dist[k]=dist[k]*((oxygen+polarconfig::mapunit[i][j+8])->charge[k])/2.0;

			}
			sum_together(sum,dist,3);
			delete [] dist;
		}
		for(size_t j=0;j<8;j++){
			dist=distance(B+i,A+polarconfig::mapunit[i][j],p);
			for(size_t k=0;k<3;k++){
				//now this guy turn into polar.
				dist[k]=dist[k]*((A+polarconfig::mapunit[i][j])->charge[k])/8.0;
			}
			sum_together(sum,dist,3);
			delete [] dist;
      }
    for(size_t m=0;m<3;m++){
      sum[m]=sum[m]/volume*16;
    }
    offset=3*sizeof(double)*i;
    MPI_File_write_at_all(fh,initial_offset+offset,sum,3,MPI_DOUBLE,&status);
    }
    MPI_Barrier(MPI_COMM_WORLD);
	}
  MPI_File_close(&fh);
}
void dispA_calculate_search(atom *A,atom *B,atom *oxygen,double *p,int Nx,int Ny,int Nz){
	double* dist;
	double* sum=new double[3];
  int world_rank;
  int world_size;
  MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
  MPI_Comm_size(MPI_COMM_WORLD,&world_size);
  int MPI_LOOP_COUNT=ceil((Nx*Ny*Nz+0.0)/world_size);
  int i=0;
  MPI_File fh;
  MPI_File_open(MPI_COMM_WORLD,"dispA.bin", MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_APPEND,MPI_INFO_NULL,&fh);
  MPI_Offset initial_offset;
  MPI_Offset offset;
  MPI_File_get_position(fh,&initial_offset);
  MPI_Status status;
  int* neighbor;
	for(size_t layer=0;layer<MPI_LOOP_COUNT;layer++){
    i=layer*world_size+world_rank;
    if(i<Nx*Ny*Nz){
		for(size_t k=0;k<3;k++){
			sum[k]=0.0;
		}
		for(size_t j=0;j<12;j++){
			dist=distance(A+i,polarconfig::mapunitA[i][j]+oxygen,p);
			sum_together(sum,dist,3);
			delete [] dist;
		}
    for(size_t m=0;m<3;m++){
      sum[m]=sum[m]/12.0;
    }
    offset=3*sizeof(double)*i;
    MPI_File_write_at_all(fh,initial_offset+offset,sum,3,MPI_DOUBLE,&status);
    }
    MPI_Barrier(MPI_COMM_WORLD);
	}
  MPI_File_close(&fh);
}
