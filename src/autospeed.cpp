#include "autospeed.h"
#include <list>
#include <iostream>
#include <fstream>
#include <complex.h>
#include <fftw3.h>
double innerprod(std::list<double*>& ve_list,int cell,size_t framei,size_t framej){
	double sum=0.0;
	std::list<double*>::iterator framei_iterator=ve_list.begin();
	std::list<double*>::iterator framej_iterator=ve_list.begin();
	for(size_t i=0;i<framei;i++){
			framei_iterator++;
	}
	for(size_t j=0;j<framej;j++){
		framej_iterator++;
	}
	for(size_t i=0;i<5*3*cell*cell*cell;i++){
		sum=sum+(*framei_iterator)[i]*(*framej_iterator)[i];
	}
	return sum;
}
double autocorre(std::list<double*>& ve_list,int cell,size_t interval){
	double sum=0.0;
	size_t len=ve_list.size();
	for(size_t i=0;i<len-interval;i++){
		sum=sum+innerprod(ve_list,cell,i,i+interval);
	}
	return sum/(len-interval);
}
void autospeed(std::list<double*>& ve_list,int cell){
	size_t len=ve_list.size();
	/*
	 * this part is used to test whether the output is what we want.
	for(std::list<double*>::iterator i=ve_list.begin();i!=ve_list.end();i++){
		for(int j=0;j<5*3*cell*cell*cell;j++){
				if(j%3==0){
					std::cout<<std::endl;
				}
				std::cout<<*((*i)+j)<<" ";
		}
	}
	*/
	std::fstream autoout;
	autoout.open("autocorrelation_of_velocity.txt",std::fstream::out);
	for(size_t i=0;i<len-1;i++){
		autoout<<autocorre(ve_list,cell,i)<<std::endl;
	}
}
