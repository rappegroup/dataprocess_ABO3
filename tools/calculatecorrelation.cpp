#include <mpi.h>
#include <string>
#include <cmath>
#include <fstream>
#include <sstream>
#include <list>
int* changeindex(int index,int cell){
	int* re=new int[3];
	re[2]=floor(index/(cell*cell));
	index=index-re[2]*cell*cell;
	re[1]=floor(index/cell);
	re[0]=index-cell*re[1];
	return re;
}
int main(int argc,char* argv[]){
	std::fstream fs;
	double celllength=4.04;/*cell length*/
	double cutoff=100;/*angstrom*/
	double distance=0.0;
	int cell=72;
	int frame=100000/200;
	std::string temp;
	std::stringstream ss;
	std::list<int> atomlist;
	int atomid;
	MPI_Init(NULL,NULL);
	int world_rank,world_size;
	MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
	MPI_Comm_size(MPI_COMM_WORLD,&world_size);
	std::list<int> pairlist;
	int* pointA;
	int* pointB;
	if(world_size==0){
	fs.open(argv[1],std::fstream::in);
	while(getline(fs,temp)){
		ss.clear();
		ss.str(temp);
		ss>>atomid;
		atomlist.push_back(atomid);
	}
	fs.close();
	int* atomref=new int [atomlist.size()];
	atomlist.sort([](int a,int b)->bool{ return a<b;});
	/*searching the pairs that satifying the cutoff criterial*/
	for(std::list<int>::iterator a=atomlist.begin();a!=atomlist.end();a++)
		for(std::list<int>::iterator b=atomlist.begin();b!=atomlist.end();b++){
			pointA=changeindex(*a,cell);
			pointB=changeindex(*b,cell);
			distance=0.0;
			for(size_t i=0;i<3;i++){
				distance=(pointA[i]-pointB[i])*(pointA[i]-pointB[i])+distance;
			}
			distance=sqrt(distance);
			if((*a < *b)&&( distance < cutoff /celllength )){
				pairlist.push_back(*a);
				pairlist.push_back(*b);
			}
			delete [] pointA;
			delete [] pointB;
		}
	std::cout<<"there are "<<pairlist.size()/2<<" pairs"<<std::endl;
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
	MPI_Bcast(pairA,pairsize,MPI::INT,0,MPI_COMM_WORLD);
	MPI_Bcast(pairB,pairsize,MPI::INT,0,MPI_COMM_WORLD);
	MPI_File fh;
	MPI_File_open(MPI_COMM_WORLD,"polar_direction_Asite.bin",MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);
	/*pair pattern (distance,timedelay correlation)*/
	double sum=0.0;
	MPI_Offset offsetone,offsettwo;
	MPI_Status status;
	double directA[3]={0.0,0.0,0.0};
	double directB[3]={0.0,0.0,0.0};
	double* timedelay=new double [frame];
	MPI_File correlation;
	MPI_File_open(MPI_COMM_WORLD,"polar_correlation.bin",MPI_MODE_WRONLY,MPI_INFO_NULL,&correlation);
	for(size_t i=world_rank;i<pairsize;i=i+world_size){
		for(size_t j=0;j<frame;j++){
			/*the time delay*/
			sum=0.0;
		  for(size_t k=0;k<frame;k++){
				/*(k,k+j)*/
				if(k+j>frame){
					sum=sum+0.0;
				}
				else{
					offsetone=(k*cell*cell*cell+pairA[i])*3;
					offsettwo=((k+j)*cell*cell*cell+pairB[i])*3;
					MPI_File_read_at(fh,offsetone,directA,3,MPI::DOUBLE,&status);
					MPI_File_read_at(fh,offsettwo,directB,3,MPI::DOUBLE,&status);
					for(size_t t=0;t<3;t++){
						sum=sum+directA[t]*directB[t];
					}
				}
			}
			sum=sum/frame;
			timedelay[j]=sum;
		}
		pointA=changeindex(pairA[i],cell);
		pointB=changeindex(pairB[i],cell);
		distance=0.0;
		for(size_t t=0;t<3;t++){
			distance=distance+(pointA[t]-pointB[t])*(pointA[t]-pointB[t]);
		}
		distance=sqrt(distance);
		distance=distance*celllength;
		/*store the data as (distance,timedelay correaltion)--->N+1*/
		offsetone=i*(1+frame);
		MPI_File_write_at_all(correlation,offsetone,&distance,1,MPI::DOUBLE,&status);
		offsetone=i*(1+frame)+1;
		MPI_File_write_at_all(correlation,offsetone,timedelay,frame,MPI::DOUBLE,&status);
	}
	MPI_File_close(&fh);
	MPI_File_close(&correlation);
	MPI_Finalize();
}
