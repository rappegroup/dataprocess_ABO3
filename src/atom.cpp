#include "atom.h"
#include "polarconfig.h"
#include <iostream>
#include <sstream>
#include <cmath>
#include <fstream>
#include <vector>
#include <new>
#include <map>
#include <list>
#include <math.h>
#include <mpi.h>
#include <algorithm>
//a=origin,b=end
double* distance(atom* a,atom* b,double* p){
	double* dist=new double[3];
	double temp;
	for(size_t i=0;i<3;i++){
		temp=b->position[i]-a->position[i];
		temp=(temp/p[i]-round(temp/p[i]))*p[i];
		dist[i]=temp;
	}
	return dist;
}

double* distance(double* a,double* b,double* p){
	double* dist=new double[3];
	double temp;
	for(size_t i=0;i<3;i++){
		temp=b[i]-a[i];
		temp=(temp/p[i]-round(temp/p[i]))*p[i];
		dist[i]=temp;
	}
	return dist;
}
/*no periodical boudary condition*/
double* distance(double* a,double* b){
    double* dist=new double[3];
    double temp;
    for(size_t i=0;i<3;i++){
        temp=a[i]-b[i];
        dist[i]=temp;
    }
    return dist;
}
//compute the distance from a to b
double far(atom* a,atom* b,double* p){
	double* temp;
	double sum=0;
	temp=distance(a,b,p);
	for(size_t i=0;i<3;i++){
		sum=sum+temp[i]*temp[i];
	}
  delete [] temp;
	return sqrt(sum);
}
double far(double* a,double* b,double* p){
  double* temp;
  double sum=0.0;
  temp=distance(a,b,p);
  for(size_t i=0;i<3;i++){
    sum=sum+temp[i]*temp[i];
  }
  delete [] temp;
  return sqrt(sum);
}
int* changeindex(int index,int cell){
	int* re=new int[3];
	re[2]=floor(index/(cell*cell));
	index=index-re[2]*cell*cell;
	re[1]=floor(index/cell);
	re[0]=index-cell*re[1];
	return re;
}
void sort(double* input,int dim){
	double good;
	for(size_t i=0;i<dim-1;i++){
		for(size_t j=0;j<dim-i-1;j++){
			if(input[j]>input[j+1]){
			good=input[j];
			input[j]=input[j+1];
			input[j+1]=good;
			}
			else continue;
		}
	}
}
int changeback(int x,int y, int z,int cell){
	return (x+cell)%cell+(y+cell)%cell*cell+(z+cell)%cell*cell*cell;
}
int* neighbor_o_forB(int index,int cell){
	int* index_3D=changeindex(index,cell);
	int* nei=new int[6];
	nei[0]=changeback(index_3D[0],index_3D[1],index_3D[2],cell);
	nei[1]=changeback(index_3D[0],index_3D[1],index_3D[2]+1,cell);
	nei[2]=changeback(index_3D[0],index_3D[1],index_3D[2],cell)+cell*cell*cell;
	nei[3]=changeback(index_3D[0],index_3D[1]+1,index_3D[2],cell)+cell*cell*cell;
	nei[4]=changeback(index_3D[0],index_3D[1],index_3D[2],cell)+2*cell*cell*cell;
	nei[5]=changeback(index_3D[0]+1,index_3D[1],index_3D[2],cell)+2*cell*cell*cell;
	return nei;
}
int* neighbor_o_forA(int index,int cell){
	int* index_3D=changeindex(index,cell);
	int* nei=new int[12];
	nei[0]=changeback(index_3D[0],index_3D[1],index_3D[2],cell);
	nei[1]=changeback(index_3D[0]-1,index_3D[1],index_3D[2],cell);
	nei[2]=changeback(index_3D[0],index_3D[1]-1,index_3D[2],cell);
	nei[3]=changeback(index_3D[0]-1,index_3D[1]-1,index_3D[2],cell);
	nei[4]=changeback(index_3D[0],index_3D[1],index_3D[2],cell)+cell*cell*cell;
	nei[5]=changeback(index_3D[0]-1,index_3D[1],index_3D[2],cell)+cell*cell*cell;
	nei[6]=changeback(index_3D[0]-1,index_3D[1],index_3D[2]-1,cell)+cell*cell*cell;
	nei[7]=changeback(index_3D[0],index_3D[1],index_3D[2]-1,cell)+cell*cell*cell;
	nei[8]=changeback(index_3D[0],index_3D[1],index_3D[2],cell)+2*cell*cell*cell;
	nei[9]=changeback(index_3D[0],index_3D[1]-1,index_3D[2],cell)+2*cell*cell*cell;
	nei[10]=changeback(index_3D[0],index_3D[1],index_3D[2]-1,cell)+2*cell*cell*cell;
	nei[11]=changeback(index_3D[0],index_3D[1]-1,index_3D[2]-1,cell)+2*cell*cell*cell;
	return nei;
}
int* neighbor_A_forB(int index,int cell){
	int* index_3D=changeindex(index,cell);
	int* nei=new int[8];
	nei[0]=changeback(index_3D[0],index_3D[1],index_3D[2],cell);
	nei[1]=changeback(index_3D[0]+1,index_3D[1],index_3D[2],cell);
	nei[2]=changeback(index_3D[0],index_3D[1]+1,index_3D[2],cell);
	nei[3]=changeback(index_3D[0]+1,index_3D[1]+1,index_3D[2],cell);
	nei[4]=changeback(index_3D[0],index_3D[1],index_3D[2]+1,cell);
	nei[5]=changeback(index_3D[0]+1,index_3D[1],index_3D[2]+1,cell);
	nei[6]=changeback(index_3D[0],index_3D[1]+1,index_3D[2]+1,cell);
	nei[7]=changeback(index_3D[0]+1,index_3D[1]+1,index_3D[2]+1,cell);
	return nei;
}
void sum_together(double* sum,double* add,int len){
	for(size_t i=0;i<len;i++){
		sum[i]=sum[i]+add[i];
	}
}
double average(std::list<double> &input){
	double sum=0.0;
	for(std::list<double>::iterator a=input.begin();a!=input.end();a++){
		sum=sum+*a;
	}
	return sum/input.size();
}
double sum_together(std::list<double> &input){
  double sum=0.0;
  for(std::list<double>::iterator a=input.begin();a!=input.end();a++){
    sum=sum+*a;
  }
  return sum;
}
double variance(std::list<double> &input){
	double sum=0.0;
	double ave=average(input);
	for(std::list<double>::iterator a=input.begin();a!=input.end();a++){
		sum=sum+(*a-ave)*(*a-ave);
	}
	return sum/input.size();
}
double* polar_average(atom *A,atom *B,atom *oxygen,double *p,int cell){
	std::list<double> px;
	std::list<double> py;
	std::list<double> pz;
	double volume=1.0;
	for(size_t k=0;k<3;k++){
		volume=volume*(p[k]/cell);
	}
	int* neighbor;
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
	for(size_t layer=0;layer<MPI_LOOP_COUNT;layer++){
    i=layer*world_size+world_rank;
    if(i<cell*cell*cell){
		neighbor=neighbor_o_forB(i,cell);
		for(size_t k=0;k<3;k++){
			sum[k]=0.0;
		}
		for(size_t j=0;j<6;j++){
			dist=distance(B+i,neighbor[j]+oxygen,p);
			for(size_t k=0;k<3;k++){
				dist[k]=dist[k]*((oxygen+neighbor[j])->charge[k])/2.0;
			}
			sum_together(sum,dist,3);
			delete [] dist;
		}
		delete [] neighbor;
		neighbor=neighbor_A_forB(i,cell);
		for(size_t j=0;j<8;j++){
			dist=distance(B+i,A+neighbor[j],p);
			for(size_t k=0;k<3;k++){
				//now this guy turn into polar.
				dist[k]=dist[k]*((A+neighbor[j])->charge[k])/8.0;
			}
			sum_together(sum,dist,3);
			delete [] dist;
      }
		px.push_back(sum[0]/volume*16);//16 is aim at converting the units from e to C
		py.push_back(sum[1]/volume*16);//16 is aim at converting the units from e to C
		pz.push_back(sum[2]/volume*16);//16 is aim at converting the units from e to C
		delete [] neighbor;
    for(size_t m=0;m<3;m++){
      sum[m]=sum[m]/volume*16;
    }
    offset=3*sizeof(double)*i;
    MPI_File_write_at_all(fh,initial_offset+offset,sum,3,MPI_DOUBLE,&status);
    }
    MPI_Barrier(MPI_COMM_WORLD);
	}
  MPI_File_close(&fh);
  sum[0]=sum_together(px);
  sum[1]=sum_together(py);
  sum[2]=sum_together(pz);
	double* reduce=new double[3];
  MPI_Reduce(sum,reduce,3,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  delete [] sum;
  for(size_t i=0;i<3;i++){
    reduce[i]=reduce[i]/cell/cell/cell;
  }
	return reduce;
}
void calculate_local_die(int cell,double local_volume,double temperature){
  int world_rank;
  int world_size;
  MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
  MPI_Comm_size(MPI_COMM_WORLD,&world_size);
  int MPI_LOOP_COUNT=ceil((cell*cell*cell+0.0)/world_size);
  MPI_File fh,fdieout;
  MPI_File_open(MPI_COMM_WORLD,"local_polar.bin",MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);
  MPI_File_open(MPI_COMM_WORLD,"local_die.bin",MPI_MODE_CREATE | MPI_MODE_WRONLY,MPI_INFO_NULL,&fdieout);
  MPI_Offset offset;
  MPI_File_get_size(fh,&offset);
  MPI_Status status;
  int simulation_time=offset/sizeof(double)/3/cell/cell/cell;
  double* px=new double[simulation_time];
  double* py=new double[simulation_time];
  double* pz=new double[simulation_time];
  double avepx,avepy,avepz;
  double* epx=new double[MPI_LOOP_COUNT];
  double* epy=new double[MPI_LOOP_COUNT];
  double* epz=new double[MPI_LOOP_COUNT];
  double diex;
  double diey;
  double diez;
  int i;
  for(size_t layer=0;layer<MPI_LOOP_COUNT;layer++){
    i=layer*world_size+world_rank;
    if(i<cell*cell*cell){
    avepx=0.0;
    avepy=0.0;
    avepz=0.0;
    for(size_t k=0;k<simulation_time;k++){
      offset=(3*(k*cell*cell*cell+i))*sizeof(double);
      MPI_File_read_at(fh,offset,&px[k],1,MPI_DOUBLE,&status);
      offset=offset+sizeof(double);
      MPI_File_read_at(fh,offset,&py[k],1,MPI_DOUBLE,&status);
      offset=offset+sizeof(double);
      MPI_File_read_at(fh,offset,&pz[k],1,MPI_DOUBLE,&status);
      avepx=avepx+px[k];
      avepy=avepy+py[k];
      avepz=avepz+pz[k];
    }
    avepx=avepx/simulation_time;
    avepy=avepy/simulation_time;
    avepz=avepz/simulation_time;
    diex=0.0;
    diey=0.0;
    diez=0.0;
    for(size_t k=0;k<simulation_time;k++){
      diex=(px[k]-avepx)*(px[k]-avepx)+diex;
      diey=(py[k]-avepy)*(py[k]-avepy)+diey;
      diez=(pz[k]-avepz)*(pz[k]-avepz)+diez;
    }
    diex=diex/simulation_time;
    diey=diey/simulation_time;
    diez=diez/simulation_time;
    diex=dielectric(diex,local_volume,temperature);
    diey=dielectric(diey,local_volume,temperature);
    diez=dielectric(diez,local_volume,temperature);
    offset=i*3*sizeof(double);
    MPI_File_write_at_all(fdieout,offset,&diex,1,MPI_DOUBLE,&status);
    offset=offset+sizeof(double);
    MPI_File_write_at_all(fdieout,offset,&diey,1,MPI_DOUBLE,&status);
    offset=offset+sizeof(double);
    MPI_File_write_at_all(fdieout,offset,&diez,1,MPI_DOUBLE,&status);
    }
  }
  delete [] px;
  delete [] py;
  delete [] pz;
  delete [] epx;
  delete [] epy;
  delete [] epz;
}
/*Asite displacement vector:*/
void displace_A_unit(atom* A,atom* oxygen,double* p,int cell){
  int* neighbor;
	double* sum=new double[3];
  double* dist;
  double sumall;
  int world_rank;
  int world_size;
  MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
  MPI_Comm_size(MPI_COMM_WORLD,&world_size);
  int MPI_LOOP_COUNT=ceil((cell*cell*cell+0.0)/world_size);
  int i=0;
  MPI_File fh;
  MPI_Offset initial_offset,offset;
  MPI_Status status;
  MPI_File_open(MPI_COMM_WORLD,"polar_direction_Asite.bin",MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_APPEND,MPI_INFO_NULL,&fh);
  MPI_File_get_size(fh,&initial_offset);
	for(size_t layer=0;layer<MPI_LOOP_COUNT;layer++){
    i=layer*world_size+world_rank;
    if(i<cell*cell*cell){
    offset=initial_offset+i*3*sizeof(double);
		neighbor=neighbor_o_forA(i,cell);
		for(size_t k=0;k<3;k++){
			sum[k]=0.0;
		}
		for(size_t t=0;t<12;t++){
			dist=distance(A+i,neighbor[t]+oxygen,p);
			sum_together(sum,dist,3);
      delete [] dist;
		} 
    delete [] neighbor;
   sumall=0.0;
   for(size_t k=0;k<3;k++){
     sumall=sum[k]*sum[k]+sumall;
   }
   for(size_t k=0;k<3;k++){
     sum[k]=sum[k]/sqrt(sumall);
   }
   MPI_File_write_at_all(fh,offset,sum,3,MPI_DOUBLE,&status);
    }
	}
  delete [] sum;
  MPI_File_close(&fh);
}
/*displace Bsite unit vector*/
void displace_B_unit(atom* B,atom* oxygen,double* p,int cell){
	int* neighbor;
	double* dist;
	double* sum=new double[3];
  int world_rank;
  int world_size;
  MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
  MPI_Comm_size(MPI_COMM_WORLD,&world_size);
  int MPI_LOOP_COUNT=ceil((cell*cell*cell+0.0)/world_size);
  int i=0;
  double sumall=0.0;
  MPI_File fh;
  MPI_Offset initial_offset,offset;
  MPI_Status status;
  MPI_File_open(MPI_COMM_WORLD,"polar_direction_Bsite.bin",MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_APPEND,MPI_INFO_NULL,&fh);
  MPI_File_get_size(fh,&initial_offset);
	for(size_t layer=0;layer<MPI_LOOP_COUNT;layer++){
    i=layer*world_size+world_rank;
    if(i<cell*cell*cell){
    offset=initial_offset+i*3*sizeof(double);
		neighbor=neighbor_o_forB(i,cell);
		for(size_t k=0;k<3;k++){
			sum[k]=0.0;
		}
		for(size_t j=0;j<6;j++){
			dist=distance(B+i,neighbor[j]+oxygen,p);
		  sum_together(sum,dist,3);
      delete [] dist;
		}
    sumall=0.0;
    for(size_t j=0;j<3;j++){
      sumall=sum[j]*sum[j]+sumall;
    }
    for(size_t j=0;j<3;j++){
      sum[j]=sum[j]/sqrt(sumall);
    }
    MPI_File_write_at_all(fh,offset,sum,3,MPI_DOUBLE,&status);
    }
	}
  delete [] sum;
  MPI_File_close(&fh);
}
/*A more general way for Displacement Calculation*/
double* displace_average_Asite(atom* A,atom* oxygen,double* p,int cell,char type_id){
  std::list<double> dx;
  std::list<double> dy;
  std::list<double> dz;
  int* neighbor;
	double* sum=new double[3];
  double* dist=new double[3];
  int world_rank;
  int world_size;
  MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
  MPI_Comm_size(MPI_COMM_WORLD,&world_size);
  int MPI_LOOP_COUNT=ceil((cell*cell*cell+0.0)/world_size);
  int i=0;
  int count_type=0;
	for(size_t layer=0;layer<MPI_LOOP_COUNT;layer++){
    i=layer*world_size+world_rank;
    if(i<cell*cell*cell){
		if(A[i].type==type_id){
    count_type=count_type+1;
		neighbor=neighbor_o_forA(i,cell);
		for(size_t k=0;k<3;k++){
			sum[k]=0.0;
		}
		for(size_t t=0;t<12;t++){
			dist=distance(A+i,neighbor[t]+oxygen,p);
			sum_together(sum,dist,3);
		}
    delete [] neighbor;
		dx.push_back(sum[0]/12.0);
		dy.push_back(sum[1]/12.0);
		dz.push_back(sum[2]/12.0);
		}
    }
	}
  int total_count=0;
  MPI_Reduce(&count_type,&total_count,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  delete [] sum;
	double* dm=new double[3];
	dm[0]=sum_together(dx);
	dm[1]=sum_together(dy);
	dm[2]=sum_together(dz);
  double* reduce=new double [3];
  MPI_Reduce(dm,reduce,3,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  delete [] dm;
  for(size_t i=0;i<3;i++){
    reduce[i]=reduce[i]/total_count;
  }
	return reduce;
}
double* displace_average_Bsite(atom* B,atom* oxygen,double* p,int cell,char type_id){
	std::list<double> dx;
	std::list<double> dy;
	std::list<double> dz;
	int* neighbor;
	double* dist;
	double* sum=new double[3];
  int world_rank;
  int world_size;
  MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
  MPI_Comm_size(MPI_COMM_WORLD,&world_size);
  int MPI_LOOP_COUNT=ceil((cell*cell*cell+0.0)/world_size);
  int i=0;
  int count_type=0;
	for(size_t layer=0;layer<MPI_LOOP_COUNT;layer++){
    i=layer*world_size+world_rank;
    if(i<cell*cell*cell){
    if(B[i].type==type_id){
    count_type=count_type+1;
		neighbor=neighbor_o_forB(i,cell);
		for(size_t k=0;k<3;k++){
			sum[k]=0.0;
		}
		for(size_t j=0;j<6;j++){
			dist=distance(B+i,neighbor[j]+oxygen,p);
		  sum_together(sum,dist,3);
		}
		dx.push_back(sum[0]/6.0);
		dy.push_back(sum[1]/6.0);
		dz.push_back(sum[2]/6.0);
    }
    }
	}
  int total_count=0;
  MPI_Reduce(&count_type,&total_count,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	double* dm=new double[3];
	dm[0]=sum_together(dx);
	dm[1]=sum_together(dy);
	dm[2]=sum_together(dz);
  double* reduce=new double [3];
  MPI_Reduce(dm,reduce,3,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  for(size_t i=0;i<3;i++){
    reduce[i]=reduce[i]/total_count;
  }
  delete [] dm;
	return reduce;
}
double displace_average_Asite_scalar(atom* A,atom* oxygen,double *p,int cell,char type_id){
	std::list<double> dall;
	int* neighbor;
	double* dist;
	double* sum=new double[3];
	double all=0;
  int world_rank;
  int world_size;
  MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
  MPI_Comm_size(MPI_COMM_WORLD,&world_size);
  int MPI_LOOP_COUNT=ceil((cell*cell*cell+0.0)/world_size);
  int i=0;
  int count_type=0;
	for(size_t layer=0;layer<MPI_LOOP_COUNT;layer++){
    i=layer*world_size+world_rank;
    if(i<cell*cell*cell){
		if(A[i].type==type_id){
    count_type=count_type+1;
		neighbor=neighbor_o_forA(i,cell);
		for(size_t k=0;k<3;k++){
			sum[k]=0.0;
		}
		for(size_t j=0;j<12;j++){
			dist=distance(A+i,neighbor[j]+oxygen,p);
			sum_together(sum,dist,3);
		}
    delete [] neighbor;
		all=0.0;
		for(size_t j=0;j<3;j++){
			all=all+sum[j]/12.0*sum[j]/12.0;
		}
		dall.push_back(sqrt(all));
		}
    }
	}
  int total_count=0;
  MPI_Reduce(&count_type,&total_count,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  double dsum=sum_together(dall);
  double reduce_dall;
  MPI_Reduce(&dsum,&reduce_dall,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  delete [] sum;
  reduce_dall=reduce_dall/total_count;
	return reduce_dall;
}
double displace_average_Bsite_scalar(atom* B,atom* oxygen,double* p,int cell,char type_id){
	std::list<double> dall;
	int* neighbor;
	double* dist;
	double all;
	double* sum=new double[3];
  int world_rank;
  int world_size;
  MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
  MPI_Comm_size(MPI_COMM_WORLD,&world_size);
  int MPI_LOOP_COUNT=ceil((cell*cell*cell+0.0)/world_size);
   int i=0;
   int count_type=0;
	for(size_t layer=0;layer<MPI_LOOP_COUNT;layer++){
    i=layer*world_size+world_rank;
    if(i<cell*cell*cell){
    if(B[i].type==type_id){
    count_type++;
		neighbor=neighbor_o_forB(i,cell);
		for(size_t k=0;k<3;k++){
			sum[k]=0.0;
		}
		for(size_t j=0;j<6;j++){
			dist=distance(B+i,neighbor[j]+oxygen,p);
		  sum_together(sum,dist,3);
		}
    delete [] neighbor;
		all=0.0;
		for(size_t k=0;k<3;k++){
			all=all+sum[k]/6.0*sum[k]/6.0;
		}
		dall.push_back(sqrt(all));
    }
    }
	}
  int total_count=0;
  MPI_Reduce(&count_type,&total_count,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  double dsum=sum_together(dall);
  double reduce_dall=0;
  MPI_Reduce(&dsum,&reduce_dall,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  delete [] sum;
  reduce_dall=reduce_dall/total_count;
	return reduce_dall;
}
/*End General way of doing this*/
//compute the angle a---b----c
double tiltangle(atom* a,atom* b,atom* c,double* p){
	double ab=far(a,b,p);
	double bc=far(b,c,p);
	double ac=far(a,c,p);
	double theta;
	theta=(ab*ab+bc*bc-ac*ac)/2/ab/bc;
	if(theta>-1 && theta< 1){
	   theta=(180-acos(theta)/(3.141592653)*180.0);
	}
	else{
		 theta=0;
	}
	return fabs(theta);
}

double displace_average_B_scalar(atom* B,atom* oxygen,double* p,int cell){
	std::list<double> dall;
	int* neighbor;
	double* dist;
	double all;
  int world_rank;
  int world_size;
  MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
  MPI_Comm_size(MPI_COMM_WORLD,&world_size);
  int MPI_LOOP_COUNT=ceil((cell*cell*cell+0.0)/world_size);
  int i=0;
	double* sum=new double[3];
	for(size_t layer=0;layer<MPI_LOOP_COUNT;layer++){
    i=layer*world_size+world_rank;
    if(i<cell*cell*cell){
		neighbor=neighbor_o_forB(i,cell);
		for(size_t k=0;k<3;k++){
			sum[k]=0.0;
		}
		for(size_t j=0;j<6;j++){
			dist=distance(B+i,neighbor[j]+oxygen,p);
		  sum_together(sum,dist,3);
			delete [] dist;
		}
		all=0.0;
		for(size_t k=0;k<3;k++){
			all=all+sum[k]/6.0*sum[k]/6.0;
		}
		delete [] neighbor;
		dall.push_back(sqrt(all));
    }
	}
  double displace_all;
  double displace_reduce;
  displace_all=sum_together(dall);
  MPI_Reduce(&displace_all,&displace_reduce,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	delete [] sum;
	return displace_reduce/cell/cell/cell;
}
double* displace_average_B(atom* B,atom* oxygen,double* p,int cell){
	std::list<double> dx;
	std::list<double> dy;
	std::list<double> dz;
	int* neighbor;
	double* dist;
	double* sum=new double[3];
  int world_rank;
  int world_size;
  MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
  MPI_Comm_size(MPI_COMM_WORLD,&world_size);
  int MPI_LOOP_COUNT=ceil((cell*cell*cell+0.0)/world_size);
	int i=0;
  for(size_t layer=0;layer<MPI_LOOP_COUNT;layer++){
    i=layer*world_size+world_rank;
    if(i<cell*cell*cell){
		neighbor=neighbor_o_forB(i,cell);
		for(size_t k=0;k<3;k++){
			sum[k]=0.0;
		}
		for(size_t j=0;j<6;j++){
			dist=distance(B+i,neighbor[j]+oxygen,p);
		  sum_together(sum,dist,3);
			delete [] dist;
		}
		dx.push_back(sum[0]/6.0);
		dy.push_back(sum[1]/6.0);
		dz.push_back(sum[2]/6.0);
		delete [] neighbor;
    }
	}
  sum[0]=sum_together(dx);
  sum[1]=sum_together(dy);
  sum[2]=sum_together(dz);
  double* reduce=new double[3];
  MPI_Reduce(sum,reduce,3,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	delete [] sum;
  for(size_t i=0;i<3;i++){
    reduce[i]=reduce[i]/cell/cell/cell;
  }
	return reduce;
}
double norm(double* p,int dim){
    double sum=0.0;
    for(size_t i=0;i<dim;i++){
        sum=sum+p[i]*p[i];
    }
    return sqrt(sum);
}
void analyzepolar(atom* A,atom* B,atom* oxygen,double* period,int cell){
		double* dispba;
	    double* dispca;
	    double* dispB;
		double* polar;
		double disp_scalar;
    	int* index;
    	int a;
    	int b;
    	int c;
    	double angle;
	    dispB=displace_average_B(B,oxygen,period,cell);
		disp_scalar=displace_average_B_scalar(B,oxygen,period,cell);
		polarconfig::disp_B_scalar.push_back(disp_scalar);
		polarconfig::disp_allB_x.push_back(dispB[0]);
		polarconfig::disp_allB_y.push_back(dispB[1]);
		polarconfig::disp_allB_z.push_back(dispB[2]);
		polar=polar_average(A,B,oxygen,period,cell);
		polarconfig::px.push_back(polar[0]);
		polarconfig::py.push_back(polar[1]);
		polarconfig::pz.push_back(polar[2]);
			//compute the tilt angle now;
      MPI_Barrier(MPI_COMM_WORLD);
      int world_rank;
      int world_size;
      MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
      MPI_Comm_size(MPI_COMM_WORLD,&world_size);
      int MPI_LOOP_COUNT=ceil((cell*cell*cell+0.0)/world_size);
      int i=0; 
       int count_type=0;
       std::list<double> tilt_angle_one;
       std::list<double> tilt_angle_two;
       std::list<double> tilt_angle_three;
	  for(size_t layer=0;layer<MPI_LOOP_COUNT;layer++){
        i=layer*world_size+world_rank;
        if(i<cell*cell*cell){
				index=changeindex(i,cell);
				a=changeback(index[0],index[1],index[2],cell);
				b=changeback(index[0],index[1],index[2]+1,cell);
				c=changeback(index[0],index[1],index[2]+2,cell);
				delete [] index;
				angle=tiltangle(a+oxygen,b+oxygen,c+oxygen,period);
				tilt_angle_one.push_back(angle);
        }
			}
      double tilt_local=sum_together(tilt_angle_one);
      double tilt_angle_one_sum;
      MPI_Reduce(&tilt_local,&tilt_angle_one_sum,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      tilt_angle_one_sum=tilt_angle_one_sum/cell/cell/cell;
			for(size_t layer=0;layer<MPI_LOOP_COUNT;layer++){
        i=layer*world_size+world_rank;
        if(i<cell*cell*cell){
				index=changeindex(i,cell);
				a=changeback(index[0],index[1],index[2],cell);
				b=changeback(index[0],index[1]+1,index[2],cell);
				c=changeback(index[0],index[1]+2,index[2],cell);
				delete [] index;
				angle=tiltangle(cell*cell*cell+a+oxygen,cell*cell*cell+b+oxygen,c+oxygen+cell*cell*cell,period);
				tilt_angle_two.push_back(angle);
        }
			}
      tilt_local=sum_together(tilt_angle_two);
      double tilt_angle_two_sum;
      MPI_Reduce(&tilt_local,&tilt_angle_two_sum,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      tilt_angle_two_sum=tilt_angle_two_sum/cell/cell/cell;
		for(size_t layer=0;layer<MPI_LOOP_COUNT;layer++){
        i=layer*world_size+world_rank;
        if(i<cell*cell*cell){
				index=changeindex(i,cell);
				a=changeback(index[0],index[1],index[2],cell);
				b=changeback(index[0]+1,index[1],index[2],cell);
				c=changeback(index[0]+2,index[1],index[2],cell);
				delete [] index;
				angle=tiltangle(2*cell*cell*cell+a+oxygen,2*cell*cell*cell+b+oxygen,c+oxygen+2*cell*cell*cell,period);
				tilt_angle_three.push_back(angle);
			}
      }
      tilt_local=sum_together(tilt_angle_three);
      double tilt_angle_three_sum;
      MPI_Reduce(&tilt_local,&tilt_angle_three_sum,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      tilt_angle_three_sum=tilt_angle_three_sum/cell/cell/cell;
      polarconfig::tilt_one_ave=tilt_angle_one_sum;
      polarconfig::tilt_two_ave=tilt_angle_two_sum;
      polarconfig::tilt_three_ave=tilt_angle_three_sum;
      polarconfig::tilt_angle_one.push_back(polarconfig::tilt_one_ave);
      polarconfig::tilt_angle_two.push_back(polarconfig::tilt_two_ave);
      polarconfig::tilt_angle_three.push_back(polarconfig::tilt_three_ave);
}
/*here polarvar are in units of (e/A^3)^2,volume are in units of A^3,temperature are in units of K*/
double dielectric(double polarvar,double volume,double temp){
	return 1e-30/(1.38*1e-23*8.85*1e-12)*polarvar*volume/temp;
	/*1e-30 is to convert the unit of A^3 to m^3
	 *1.38*1e-23 is kb boltzmann constant.
	 * */
}
void outpolar(){
    std::fstream fileout;
	double lx,ly,lz;
	lx=average(polarconfig::la_x);
	ly=average(polarconfig::la_y);
	lz=average(polarconfig::la_z);
	fileout.open("result.txt",std::fstream::out);
	fileout<<"the average lattice constant is:"<<std::endl;
	fileout<<lx<<" "<<ly<<" "<<lz<<std::endl;
	fileout<<"the average B site cations displacement is:"<<std::endl;
	fileout<<fabs(average(polarconfig::disp_allB_x))<<" "<<fabs(average(polarconfig::disp_allB_y))<<" "<<fabs(average(polarconfig::disp_allB_z))<<std::endl;
	fileout<<"the averaget O6 tilt angle in three direction is:"<<std::endl;
	fileout<<fabs(polarconfig::tilt_one_ave)<<" "<<fabs(polarconfig::tilt_two_ave)<<" "<<fabs(polarconfig::tilt_two_ave)<<std::endl;
	fileout<<"the scalar average B site displacement is:"<<std::endl;
	fileout<<average(polarconfig::disp_B_scalar)<<std::endl;
	std::vector<double> pall(3,0.0);
	std::vector<double> var(3,0.0);
	pall[0]=std::fabs(average(polarconfig::px));
	pall[1]=std::fabs(average(polarconfig::py));
	pall[2]=std::fabs(average(polarconfig::pz));
	std::fstream pout;
	pout.open("polar.txt",std::fstream::out);
	std::list<double>::iterator pyi=polarconfig::py.begin();
	std::list<double>::iterator pzi=polarconfig::pz.begin();
	for(std::list<double>::iterator pxi=polarconfig::px.begin();pxi!=polarconfig::px.end();pxi++){
		pout<<*(pxi)<<" "<<*pyi<<" "<<*pzi<<std::endl;
		pyi++;
		pzi++;
	}
	var[0]=variance(polarconfig::px);
	var[1]=variance(polarconfig::py);
	var[2]=variance(polarconfig::pz);
	std::map <double,double> good;
	for(size_t i=0;i<3;i++){
		good.insert(good.end(),std::pair <double,double> (pall[i],var[i]));
	}
//	sort(pall.begin(),pall.end());
	fileout<<"the polarization is (absolute value):"<<std::endl;
	fileout<<pall[0]<<" "<<pall[1]<<" "<<pall[2]<<std::endl;
	fileout<<"polarization variance is:"<<std::endl;
	fileout<<good[pall[0]]<<" "<<good[pall[1]]<<" "<<good[pall[2]]<<std::endl;
	fileout<<"the relative dielectric constant(epsilon0) is :"<<std::endl;
  int Nx=polarconfig::Nx;
  int Ny=polarconfig::Ny;
  int Nz=polarconfig::Nz;
	double dx=dielectric(good[pall[0]],lx*ly*lz*Nx*Ny*Nz,polarconfig::temperature);
	double dy=dielectric(good[pall[1]],lx*ly*lz*Nx*Ny*Nz,polarconfig::temperature);
	double dz=dielectric(good[pall[2]],lx*ly*lz*Nx*Ny*Nz,polarconfig::temperature);
	fileout<<dx<<" "<<dy<<" "<<dz<<std::endl;
	fileout<<(dx+dy+dz)/3.0<<std::endl;
	fileout.close();
    fileout.close();
}
void analyzeposition_variance(atom* A,atom* B,atom* oxygen,double* period,int cell,size_t signal){
  /*store the first position position on the disk and doing nothing*/
  int world_rank;
  int world_size;
  MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
  MPI_Comm_size(MPI_COMM_WORLD,&world_size);
  int MPI_LOOP_COUNT=ceil((cell*cell*cell+0.0)/world_size);
  MPI_Status status;
  MPI_File fh,fhA,fhB,fhAinitial,fhBinitial;
  MPI_Offset initial_offset,offset,initial_offsetA,initial_offsetB;
  double initial_position_A[3]={0.0,0.0,0.0};
  double initial_position_B[3]={0.0,0.0,0.0};
  double* bias_position;
  int i;
  if(signal==1){
    MPI_File_open(MPI_COMM_WORLD,"starting_position.A.bin",MPI_MODE_CREATE | MPI_MODE_WRONLY,MPI_INFO_NULL,&fh);
    for(size_t layer=0;layer<MPI_LOOP_COUNT;layer++){
      i=layer*world_size+world_rank;
      if(i<cell*cell*cell){
      initial_offset=i*3*sizeof(double);
      MPI_File_write_at_all(fh,initial_offset,(A+i)->position,3,MPI_DOUBLE,&status);
      }
    }
    MPI_File_close(&fh);
    MPI_File_open(MPI_COMM_WORLD,"starting_position.B.bin",MPI_MODE_CREATE | MPI_MODE_WRONLY,MPI_INFO_NULL,&fh);
    for(size_t layer=0;layer<MPI_LOOP_COUNT;layer++){
      i=layer*world_size+world_rank;
      if(i<cell*cell*cell){
      initial_offset=i*3*sizeof(double);
      MPI_File_write_at_all(fh,initial_offset,(B+i)->position,3,MPI_DOUBLE,&status);
      }
    }
    MPI_File_close(&fh);
  }
  else{
    MPI_File_open(MPI_COMM_WORLD,"starting_position.A.bin",MPI_MODE_RDONLY,MPI_INFO_NULL,&fhAinitial);
    MPI_File_open(MPI_COMM_WORLD,"starting_position.B.bin",MPI_MODE_RDONLY,MPI_INFO_NULL,&fhBinitial);
    MPI_File_open(MPI_COMM_WORLD,"position_A.bin", MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_APPEND,MPI_INFO_NULL,&fhA);
    MPI_File_open(MPI_COMM_WORLD,"position_B.bin", MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_APPEND,MPI_INFO_NULL,&fhB);
    MPI_File_get_position(fhA,&initial_offsetA);
    MPI_File_get_position(fhB,&initial_offsetB);
    for(size_t layer=0;layer<MPI_LOOP_COUNT;layer++){
    i=layer*world_size+world_rank;
    if(i<cell*cell*cell){
      offset=i*3*sizeof(double);
      MPI_File_read_at(fhAinitial,offset,initial_position_A,3,MPI_DOUBLE,&status);
      MPI_File_read_at(fhBinitial,offset,initial_position_B,3,MPI_DOUBLE,&status);
        bias_position=distance(initial_position_A,(A+i)->position,period);
      MPI_File_write_at_all(fhA,initial_offsetA+offset,bias_position,3,MPI_DOUBLE,&status);
      delete [] bias_position;
        bias_position=distance(initial_position_B,(B+i)->position,period);
      MPI_File_write_at_all(fhB,initial_offsetB+offset,bias_position,3,MPI_DOUBLE,&status);
      delete [] bias_position;
    }
    }
  MPI_File_close(&fhAinitial);
  MPI_File_close(&fhBinitial);
  MPI_File_close(&fhA);
  MPI_File_close(&fhB);

  }
 }
void calculate_local_variance(int cell,double temperature){
  int world_rank;
  int world_size;
  MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
  MPI_Comm_size(MPI_COMM_WORLD,&world_size);
  int MPI_LOOP_COUNT=ceil((cell*cell*cell+0.0)/world_size);
  MPI_Status status;
  MPI_File fh,fdieout;
  MPI_Offset offset;
  MPI_File_open(MPI_COMM_WORLD,"position_A.bin", MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);
  MPI_File_open(MPI_COMM_WORLD,"local_variance_A.bin",MPI_MODE_CREATE | MPI_MODE_WRONLY,MPI_INFO_NULL,&fdieout);
  MPI_File_get_size(fh,&offset);
  int simulation_time=offset/sizeof(double)/3/(cell*cell*cell);
  double* px=new double [simulation_time];
  double* py=new double [simulation_time];
  double* pz=new double [simulation_time];
  double avepx,avepy,avepz;
  double* epx=new double [MPI_LOOP_COUNT];
  double* epy=new double [MPI_LOOP_COUNT];
  double* epz=new double [MPI_LOOP_COUNT];
  double diex;
  double diey;
  double diez;
  int i;
  for(size_t layer=0;layer< MPI_LOOP_COUNT;layer++){
    i=layer*world_size+world_rank;
    if(i<cell*cell*cell){
    avepx=0.0;
    avepy=0.0;
    avepz=0.0;
    for(size_t k=0;k<simulation_time;k++){
      offset=(3*(k*cell*cell*cell+i))*sizeof(double);
      MPI_File_read_at(fh,offset,&px[k],1,MPI_DOUBLE,&status);
      offset=offset+sizeof(double);
      MPI_File_read_at(fh,offset,&py[k],1,MPI_DOUBLE,&status);
      offset=offset+sizeof(double);
      MPI_File_read_at(fh,offset,&pz[k],1,MPI_DOUBLE,&status);
      avepx=avepx+px[k];
      avepy=avepy+py[k];
      avepz=avepz+pz[k];
    }
    avepx=avepx/simulation_time;
    avepy=avepy/simulation_time;
    avepz=avepz/simulation_time;
    diex=0.0;
    diey=0.0;
    diez=0.0;
    for(size_t k=0;k<simulation_time;k++){
      diex=(px[k]-avepx)*(px[k]-avepx)+diex;
      diey=(py[k]-avepy)*(py[k]-avepy)+diey;
      diez=(pz[k]-avepz)*(pz[k]-avepz)+diez;
    }
    diex=diex/simulation_time/temperature;
    diey=diey/simulation_time/temperature;
    diez=diez/simulation_time/temperature;
    offset=i*3*sizeof(double);
    MPI_File_write_at_all(fdieout,offset,&diex,1,MPI_DOUBLE,&status);
    offset=offset+sizeof(double);
    MPI_File_write_at_all(fdieout,offset,&diey,1,MPI_DOUBLE,&status);
    offset=offset+sizeof(double);
    MPI_File_write_at_all(fdieout,offset,&diez,1,MPI_DOUBLE,&status);
    }
  }
  MPI_File_close(&fdieout);
  MPI_File_close(&fh);
  MPI_File_open(MPI_COMM_WORLD,"position_B.bin", MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);
  MPI_File_open(MPI_COMM_WORLD,"local_variance_B.bin",MPI_MODE_CREATE | MPI_MODE_WRONLY,MPI_INFO_NULL,&fdieout);
  MPI_File_get_size(fh,&offset);
  for(size_t layer=0;layer< MPI_LOOP_COUNT;layer++){
    i=layer*world_size+world_rank;
    if(i<cell*cell*cell){
    avepx=0.0;
    avepy=0.0;
    avepz=0.0;
    for(size_t k=0;k<simulation_time;k++){
      offset=(3*(k*cell*cell*cell+i))*sizeof(double);
      MPI_File_read_at(fh,offset,&px[k],1,MPI_DOUBLE,&status);
      offset=offset+sizeof(double);
      MPI_File_read_at(fh,offset,&py[k],1,MPI_DOUBLE,&status);
      offset=offset+sizeof(double);
      MPI_File_read_at(fh,offset,&pz[k],1,MPI_DOUBLE,&status);
      avepx=avepx+px[k];
      avepy=avepy+py[k];
      avepz=avepz+pz[k];
    }
    avepx=avepx/simulation_time;
    avepy=avepy/simulation_time;
    avepz=avepz/simulation_time;
    diex=0.0;
    diey=0.0;
    diez=0.0;
    for(size_t k=0;k<simulation_time;k++){
      diex=(px[k]-avepx)*(px[k]-avepx)+diex;
      diey=(py[k]-avepy)*(py[k]-avepy)+diey;
      diez=(pz[k]-avepz)*(pz[k]-avepz)+diez;
    }
    diex=diex/simulation_time/temperature;
    diey=diey/simulation_time/temperature;
    diez=diez/simulation_time/temperature;
    offset=i*3*sizeof(double);
    MPI_File_write_at_all(fdieout,offset,&diex,1,MPI_DOUBLE,&status);
    offset=offset+sizeof(double);
    MPI_File_write_at_all(fdieout,offset,&diey,1,MPI_DOUBLE,&status);
    offset=offset+sizeof(double);
    MPI_File_write_at_all(fdieout,offset,&diez,1,MPI_DOUBLE,&status);
    }
  }
  MPI_File_close(&fdieout);
  MPI_File_close(&fh);
}
