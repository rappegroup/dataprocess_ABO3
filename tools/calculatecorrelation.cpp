#include <mpi.h>
#include <string>
#include <cmath>
#include <fstream>
#include <sstream>
#include <list>
#include <cstdlib>
#include <iostream>
#include <time.h>
#include <stdio.h>
int* changeindex(int index,int cell){
	int* re=new int[3];
	re[2]=floor(index/(cell*cell));
	index=index-re[2]*cell*cell;
	re[1]=floor(index/cell);
	re[0]=index-cell*re[1];
	return re;
}
int main(int argc,char* argv[]){
  MPI_Init(NULL,NULL);
	int world_rank,world_size;
	MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
	MPI_Comm_size(MPI_COMM_WORLD,&world_size);
	clock_t start=clock();
	std::fstream fs;
	double celllength=4.04;/*cell length*/
	double cutin=std::atof(argv[2]);
	double cutoff=std::atof(argv[3]);/*angstrom*/
	double distance=0.0;
	double temp_double=0.0;
	int cell=48;
	int inteval=10;
	int frame=100000/200/inteval;
	std::string temp;
	std::stringstream ss;
	std::list<int> atomlist;
	int atomid;
	std::list<int> pairlist;
	int* pointA;
	int* pointB;
	if(world_rank==0){
	fs.open(argv[1],std::fstream::in);
	while(getline(fs,temp)){
		ss.clear();
		ss.str(temp);
		ss>>atomid;
		atomlist.push_back(atomid);
	}
	fs.close();
	atomlist.sort([](int a,int b)->bool{ return a<b;});
	/*searching the pairs that satifying the cutoff criterial*/
	for(std::list<int>::iterator a=atomlist.begin();a!=atomlist.end();a++)
		for(std::list<int>::iterator b=atomlist.begin();b!=atomlist.end();b++){
			pointA=changeindex(*a,cell);
			pointB=changeindex(*b,cell);
			distance=0.0;
			for(size_t i=0;i<3;i++){
				temp_double=pointA[i]-pointB[i]-round((pointA[i]-pointB[i]+0.0)/cell)*cell;
				distance=temp_double*temp_double+distance;
			}
			distance=sqrt(distance);
			if((*a < *b)&&( distance < cutoff /celllength )&&( distance > cutin /celllength )){
				pairlist.push_back(*a);
				pairlist.push_back(*b);
			}
			delete [] pointA;
			delete [] pointB;
		}
	}
	else{
	/*doing nothing and waiting for instructions*/
	}
	MPI_Barrier(MPI_COMM_WORLD);
	int pairsize;
	pairsize=pairlist.size()/2;
	MPI_Bcast(&pairsize,1,MPI::INT,0,MPI_COMM_WORLD);
	int* pairA=new int [pairsize];
	int* pairB=new int [pairsize];
	int count=0;
	for(std::list<int>::iterator a=pairlist.begin();a!=pairlist.end();a++){
		if(count%2==0){
			pairA[count/2]=*a;
		}
		else{
			pairB[(count-1)/2]=*a;
	  }
		count=count+1;
	}
	if(world_rank==0){
		std::cout<<"I have "<<pairsize<<" pairs"<<std::endl;
	std::cout<<((double)(clock()-start))/CLOCKS_PER_SEC<<" seconds"<<std::endl;
	}
	MPI_Bcast(pairA,pairsize,MPI::INT,0,MPI_COMM_WORLD);
	MPI_Bcast(pairB,pairsize,MPI::INT,0,MPI_COMM_WORLD);
	MPI_File fh;
	MPI_File_open(MPI_COMM_WORLD,"polar_direction_Asite.bin",MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);
	/*pair pattern (distance,timedelay correlation)*/
	MPI_Offset offsetone,offsettwo;
	MPI_Status status;
	double directA[3]={0.0,0.0,0.0};
	double directB[3]={0.0,0.0,0.0};
	double* timedelay=new double [frame];
	MPI_File correlation;
	MPI_File_open(MPI_COMM_WORLD,"polar_correlation.bin",MPI_MODE_CREATE | MPI_MODE_WRONLY,MPI_INFO_NULL,&correlation);
	double* correlation_function=new double [frame*pairsize];
	double* frameA=new double [cell*cell*cell*3];
	double* frameB=new double [cell*cell*cell*3];
	double* sum=new double [pairsize];
	/*parallel over the time delay*/
	for(size_t j=world_rank;j<frame;j=j+world_size){
		for(size_t k=0;k<pairsize;k++){
			sum[k]=0.0;
		}
		for(size_t k=0;k<frame;k++){
			/*(k,k+j)*/
			if(k+j>=frame){
			}
			else{
				offsetone=(k*inteval*cell*cell*cell)*3*sizeof(double);
				offsettwo=((k+j)*inteval*cell*cell*cell)*3*sizeof(double);
				MPI_File_read_at(fh,offsetone,frameA,3*cell*cell*cell,MPI::DOUBLE,&status);
				MPI_File_read_at(fh,offsettwo,frameB,3*cell*cell*cell,MPI::DOUBLE,&status);
				for(size_t m=0;m<pairsize;m++){
					for(size_t n=0;n<3;n++){
						sum[m]=sum[m]+frameA[pairA[m]*3+n]*frameB[pairB[m]*3+n];
					}
				}
			}
		}
			for(size_t m=0;m<pairsize;m++){
				sum[m]=sum[m]/frame;
			}
			offsetone=(j*pairsize)*sizeof(double);
			MPI_File_write_at(correlation,offsetone,sum,pairsize,MPI::DOUBLE,&status);
	}
	MPI_File correlation_distance;
	MPI_File_open(MPI_COMM_WORLD,"correlation_distance.bin",MPI_MODE_CREATE | MPI_MODE_WRONLY,MPI_INFO_NULL,&correlation_distance);
	for(size_t j=world_rank;j<pairsize;j=j+world_size){
		distance=0.0;
		pointA=changeindex(pairA[j],cell);
		pointB=changeindex(pairB[j],cell);
		for(size_t t=0;t<3;t++){
			temp_double=pointA[t]-pointB[t]-round((pointA[t]-pointB[t]+0.0)/cell)*cell;
			distance=distance+temp_double*temp_double;
		}
		distance=sqrt(distance);
		distance=distance*celllength;
		offsetone=j*sizeof(double);
		MPI_File_write_at(correlation_distance,offsetone,&distance,1,MPI::DOUBLE,&status);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_File_close(&fh);
	MPI_File_close(&correlation);
	MPI_File_close(&correlation_distance);
	clock_t duration=clock()-start;
	if(world_rank==0){
	std::cout<<((double)duration)/CLOCKS_PER_SEC<<" seconds"<<std::endl;
	}
	MPI_Finalize();
	std::cout<<"I finished"<<std::endl;
	exit(0);
}
