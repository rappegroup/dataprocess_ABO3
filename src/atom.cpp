#include "atom.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <map>
#include <list>
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
	for(size_t i=0;i<cell*cell*cell;i++){
		neighbor=neighbor_o_forB(i,cell);
		for(size_t k=0;k<3;k++){
			sum[k]=0.0;
		}
		for(size_t j=0;j<6;j++){
			dist=distance(B+i,neighbor[j]+oxygen,p);
			for(size_t k=0;k<3;k++){
				dist[k]=dist[k]*(oxygen->charge)/2.0;
			}
			sum_together(sum,dist,3);
		}
		neighbor=neighbor_A_forB(i,cell);
		for(size_t j=0;j<8;j++){
			dist=distance(B+i,A+neighbor[j],p);
			for(size_t k=0;k<3;k++){
				//now this guy turn into polar.
				dist[k]=dist[k]*(A->charge)/8.0;
			}
			sum_together(sum,dist,3);
		}
		px.push_back(sum[0]/volume*16);//16 is aim at converting the units from e to C
		py.push_back(sum[1]/volume*16);//16 is aim at converting the units from e to C
		pz.push_back(sum[2]/volume*16);//16 is aim at converting the units from e to C
	}
	double* pall=new double[3];
	pall[0]=average(px);
	pall[1]=average(py);
	pall[2]=average(pz);
	return pall;
}
double* displace_average_A(atom* A,atom* oxygen,double *p,int cell){
	std::list<double> dx;
	std::list<double> dy;
	std::list<double> dz;
	int* neighbor;
	double* dist;
	double* sum=new double[3];
	for(size_t i=0;i<cell*cell*cell;i++){
		neighbor=neighbor_o_forA(i,cell);
		for(size_t k=0;k<3;k++){
			sum[k]=0.0;
		}
		for(size_t j=0;j<12;j++){
			dist=distance(A+i,neighbor[j]+oxygen,p);
			sum_together(sum,dist,3);
		}
		dx.push_back(sum[0]/12.0);
		dy.push_back(sum[1]/12.0);
		dz.push_back(sum[2]/12.0);
	}
	double* dm=new double[3];
	dm[0]=average(dx);
	dm[1]=average(dy);
	dm[2]=average(dz);
	return dm;
}
double* displace_average_Ca(atom* A,atom* oxygen,double *p,int cell){
	std::list<double> dx;
	std::list<double> dy;
	std::list<double> dz;
	int* neighbor;
	double* dist;
	double* sum=new double[3];
	for(size_t i=0;i<cell*cell*cell;i++){
		if(A[i].type=='c'){
		neighbor=neighbor_o_forA(i,cell);
		for(size_t k=0;k<3;k++){
			sum[k]=0.0;
		}
		for(size_t j=0;j<12;j++){
			dist=distance(A+i,neighbor[j]+oxygen,p);
			sum_together(sum,dist,3);
		}
		dx.push_back(sum[0]/12.0);
		dy.push_back(sum[1]/12.0);
		dz.push_back(sum[2]/12.0);
		}
	}
	double* dm=new double[3];
	dm[0]=average(dx);
	dm[1]=average(dy);
	dm[2]=average(dz);
	return dm;
}
double* displace_average_Ba(atom* A,atom* oxygen,double *p,int cell){
	std::list<double> dx;
	std::list<double> dy;
	std::list<double> dz;
	int* neighbor;
	double* dist;
	double* sum=new double[3];
	for(size_t i=0;i<cell*cell*cell;i++){
		if(A[i].type=='b'){
		neighbor=neighbor_o_forA(i,cell);
		for(size_t k=0;k<3;k++){
			sum[k]=0.0;
		}
		for(size_t t=0;t<12;t++){
			dist=distance(A+i,neighbor[t]+oxygen,p);
			sum_together(sum,dist,3);
		}
		dx.push_back(sum[0]/12.0);
		dy.push_back(sum[1]/12.0);
		dz.push_back(sum[2]/12.0);
		}
	}
	double* dm=new double[3];
	dm[0]=average(dx);
	dm[1]=average(dy);
	dm[2]=average(dz);
	return dm;
}
//compute the angle a---b----c
double tiltangle(atom* a,atom* b,atom* c,double* p){
	double ab=far(a,b,p);
	double bc=far(b,c,p);
	double ac=far(a,c,p);
	double theta;
	theta=(ab*ab+bc*bc-ac*ac)/2/ab/bc;
	if(theta>-1 && theta<1){
	   theta=(180-acos(theta)/(3.141592653)*180.0);
	}
	else{
		 theta=0;
	}
	return fabs(theta);
}
//giving the scalar average version of this code
double displace_average_Ca_scalar(atom* A,atom* oxygen,double *p,int cell){
	std::list<double> dall;
	int* neighbor;
	double* dist;
	double* sum=new double[3];
	double all=0;
	for(size_t i=0;i<cell*cell*cell;i++){
		if(A[i].type=='c'){
		neighbor=neighbor_o_forA(i,cell);
		for(size_t k=0;k<3;k++){
			sum[k]=0.0;
		}
		for(size_t j=0;j<12;j++){
			dist=distance(A+i,neighbor[j]+oxygen,p);
			sum_together(sum,dist,3);
		}
		all=0.0;
		for(size_t j=0;j<3;j++){
			all=all+sum[j]/12.0*sum[j]/12.0;
		}
		dall.push_back(sqrt(all));
		}
	}
	return average(dall);
}
double displace_average_Ba_scalar(atom* A,atom* oxygen,double *p,int cell){
	std::list<double> dall;
	int* neighbor;
	double* dist;
	double* sum=new double[3];
	double all=0;
	for(size_t i=0;i<cell*cell*cell;i++){
		if(A[i].type=='b'){
		neighbor=neighbor_o_forA(i,cell);
		for(size_t k=0;k<3;k++){
			sum[k]=0.0;
		}
		for(size_t j=0;j<12;j++){
			dist=distance(A+i,neighbor[j]+oxygen,p);
			sum_together(sum,dist,3);
		}
		all=0.0;
		for(size_t j=0;j<3;j++){
			all=all+sum[j]/12.0*sum[j]/12.0;
		}
		dall.push_back(sqrt(all));
		}
	}
	return average(dall);
}
double displace_average_B_scalar(atom* B,atom* oxygen,double* p,int cell){
	std::list<double> dall;
	int* neighbor;
	double* dist;
	double all;
	double* sum=new double[3];
	for(size_t i=0;i<cell*cell*cell;i++){
		neighbor=neighbor_o_forB(i,cell);
		for(size_t k=0;k<3;k++){
			sum[k]=0.0;
		}
		for(size_t j=0;j<6;j++){
			dist=distance(B+i,neighbor[j]+oxygen,p);
		  sum_together(sum,dist,3);
		}
		all=0.0;
		for(size_t k=0;k<3;k++){
			all=all+sum[k]/6.0*sum[k]/6.0;
		}
		dall.push_back(sqrt(all));
	}
	return average(dall);
}
double* displace_average_B(atom* B,atom* oxygen,double* p,int cell){
	std::list<double> dx;
	std::list<double> dy;
	std::list<double> dz;
	int* neighbor;
	double* dist;
	double* sum=new double[3];
	for(size_t i=0;i<cell*cell*cell;i++){
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
	double* dm=new double[3];
	dm[0]=average(dx);
	dm[1]=average(dy);
	dm[2]=average(dz);
	return dm;
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
			dispba=displace_average_Ba(A,oxygen,period,cell);
			disp_scalar=displace_average_Ba_scalar(A,oxygen,period,cell);
			polarconfig::disp_ba_scalar.push_back(disp_scalar);
			polarconfig::disp_allba_x.push_back(dispba[0]);
			polarconfig::disp_allba_y.push_back(dispba[1]);
			polarconfig::disp_allba_z.push_back(dispba[2]);
			dispca=displace_average_Ca(A,oxygen,period,cell);
			disp_scalar=displace_average_Ca_scalar(A,oxygen,period,cell);
			polarconfig::disp_ca_scalar.push_back(disp_scalar);
			polarconfig::disp_allca_x.push_back(dispca[0]);
			polarconfig::disp_allca_y.push_back(dispca[1]);
			polarconfig::disp_allca_z.push_back(dispca[2]);
			polar=polar_average(A,B,oxygen,period,cell);
			//sort(polar,3);
			polarconfig::px.push_back(polar[0]);
			polarconfig::py.push_back(polar[1]);
			polarconfig::pz.push_back(polar[2]);
			//compute the tilt angle now;
			for(size_t i=0;i<cell*cell*cell;i++){
				index=changeindex(i,cell);
				a=changeback(index[0],index[1],index[2],cell);
				b=changeback(index[0],index[1],index[2]+1,cell);
				c=changeback(index[0],index[1],index[2]+2,cell);
				angle=tiltangle(a+oxygen,b+oxygen,c+oxygen,period);
				polarconfig::tilt_angle.push_back(angle);
				polarconfig::tilt_angle_one.push_back(angle);
			}
			for(size_t i=0;i<cell*cell*cell;i++){
				index=changeindex(i,cell);
				a=changeback(index[0],index[1],index[2],cell);
				b=changeback(index[0],index[1]+1,index[2],cell);
				c=changeback(index[0],index[1]+2,index[2],cell);
				angle=tiltangle(cell*cell*cell+a+oxygen,cell*cell*cell+b+oxygen,c+oxygen+cell*cell*cell,period);
				polarconfig::tilt_angle.push_back(angle);
				polarconfig::tilt_angle_two.push_back(angle);
			}
		for(size_t i=0;i<cell*cell*cell;i++){
				index=changeindex(i,cell);
				a=changeback(index[0],index[1],index[2],cell);
				b=changeback(index[0]+1,index[1],index[2],cell);
				c=changeback(index[0]+2,index[1],index[2],cell);
				angle=tiltangle(2*cell*cell*cell+a+oxygen,2*cell*cell*cell+b+oxygen,c+oxygen+2*cell*cell*cell,period);
				polarconfig::tilt_angle.push_back(angle);
				polarconfig::tilt_angle_three.push_back(angle);
			}
}
void outpolar(){
  using namespace polarconfig;
  std::fstream fileout;
	fileout.open("result.txt",std::fstream::out);
	fileout<<"the average lattice constant is:"<<std::endl;
	fileout<<average(la_x)<<" "<<average(la_y)<<" "<<average(la_z)<<std::endl;
	fileout<<"the average B site cations displacement is:"<<std::endl;
	fileout<<fabs(average(disp_allB_x))<<" "<<fabs(average(disp_allB_y))<<" "<<fabs(average(disp_allB_z))<<std::endl;
	fileout<<"the average Asite one displacement is:"<<std::endl;
	fileout<<fabs(average(disp_allba_x))<<" "<<fabs(average(disp_allba_y))<<" "<<fabs(average(disp_allba_z))<<std::endl;
	fileout<<"the average Asite two displacement is:"<<std::endl;
	fileout<<fabs(average(disp_allca_x))<<" "<<fabs(average(disp_allca_y))<<" "<<fabs(average(disp_allca_z))<<std::endl;
	fileout<<"the average O6 tilt angle is:"<<std::endl;
	fileout<<(fabs(average(tilt_angle)))<<std::endl;
	fileout<<"the averaget O6 tilt angle in three direction is:"<<std::endl;
	fileout<<(fabs(average(tilt_angle_one)))<<" "<<fabs(average(tilt_angle_two))<<" "<<fabs(average(tilt_angle_three))<<std::endl;
	fileout<<"the scalar average B site displacement is:"<<std::endl;
	fileout<<average(disp_B_scalar)<<std::endl;
	fileout<<"the scalar average Asite one dispalcement is:"<<std::endl;
	fileout<<average(disp_ba_scalar)<<std::endl;
	fileout<<"the scalar average Asite two displacement is:"<<std::endl;
	fileout<<average(disp_ca_scalar)<<std::endl;
	std::vector<double> pall(3,0.0);
	std::vector<double> var(3,0.0);
	pall[0]=std::fabs(average(px));
	pall[1]=std::fabs(average(py));
	pall[2]=std::fabs(average(pz));
	std::fstream pout;
	pout.open("polar.txt",std::fstream::out);
	std::list<double>::iterator pyi=py.begin();
	std::list<double>::iterator pzi=pz.begin();
	for(std::list<double>::iterator pxi=px.begin();pxi!=px.end();pxi++){
		pout<<*(pxi)<<" "<<*pyi<<" "<<*pzi<<std::endl;
		pyi++;
		pzi++;
	}
	var[0]=variance(px);
	var[1]=variance(py);
	var[2]=variance(pz);
	std::map <double,double> good;
	for(size_t i=0;i<3;i++){
		good.insert(good.end(),std::pair <double,double> (pall[i],var[i]));
	}
//	sort(pall.begin(),pall.end());
	fileout<<"the polarization is (absolute value):"<<std::endl;
	fileout<<pall[0]<<" "<<pall[1]<<" "<<pall[2]<<std::endl;
	fileout<<"polarization variance is:"<<std::endl;
	fileout<<good[pall[0]]<<" "<<good[pall[1]]<<" "<<good[pall[2]]<<std::endl;
	fileout.close();
}
