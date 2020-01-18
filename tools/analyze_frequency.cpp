#include <iostream>
#include <sstream>
#include <string>
#include <list>
#include <fstream>
#include <new>
#include <fftw3.h>
#include <complex>
#define PI 3.141592653
std::complex<double> dielectric_second(std::complex<double> polarvar,double volume,double temp,double frequency){
  return 1e-30/(1.38*1e-23*8.85*1e-12)*polarvar*volume/temp*sqrt(std::complex<double>(-1,0))*2.0*PI*frequency;
    /*1e-30 is to convert the unit of A^3 to m^3
     *   *1.38*1e-23 is kb boltzmann constant.
     *     * */
}
double dielectric_first(double polarvar,double volume,double temp){
   return 1e-30/(1.38*1e-23*8.85*1e-12)*polarvar*volume/temp;
    /*1e-30 is to convert the unit of A^3 to m^3
     *   *1.38*1e-23 is kb boltzmann constant.
     *     * */
}
double* polarcorrelation(double* pseries,int length){
  int outlength=length-1;
  double* correlation=new double[outlength];
  double sum=0.0;
  for(size_t i=0;i<outlength;i++){
    sum=0.0;
    for(size_t k=0;k+i<length;k++){
    sum=sum+pseries[k]*pseries[k+i];
    }
    correlation[i]=sum/(length-i);
  }
  return correlation;
}
double* polarcorrelation_zero_padding(double* pseries,int length){
  int outlength=length;
  double* correlation=new double[outlength];
  double sum=0.0;
  for(size_t i=0;i<length;i++){
    sum=0.0;
    for(size_t k=0;k+i<length;k++){
    sum=sum+pseries[k]*pseries[k+i];
    }
    correlation[i]=sum/(length);
  }
  return correlation;
}

double polarvar(double* pseries,int length){
  double sum=0.0;
  for(size_t i=0;i<length;i++){
    sum=sum+pseries[i]*pseries[i];
  }
  return sum/length;
}
double averagelist(std::list<double>& plist){
  double sum=0.0;
  for(std::list<double>::iterator a=plist.begin();a!=plist.end();a++){
    sum=sum+*a;
  }
  return sum/plist.size();
}
int main(){
  std::fstream fs;
  fs.open("polar.txt",std::fstream::in);
  std::string templine;
  std::stringstream linestream;
  double px_temp,py_temp,pz_temp;
  std::list<double> px_list,py_list,pz_list;
  int simulation_time_steps=0;
  while(getline(fs,templine)){
    linestream.clear();
    linestream.str(templine);
    linestream>>px_temp;
    linestream>>py_temp;
    linestream>>pz_temp;
    px_list.push_back(px_temp);
    py_list.push_back(py_temp);
    pz_list.push_back(pz_temp);
    simulation_time_steps++;
  }
  int equilibrium_time_steps=1000000;
  int dump_inteval=200;
  simulation_time_steps=simulation_time_steps*dump_inteval;
  int useful=(simulation_time_steps-equilibrium_time_steps)/dump_inteval;
  double* px_vector=new double[useful];
  double* py_vector=new double[useful];
  double* pz_vector=new double[useful];
  double volume=(20*4.04)*(20*4.04)*(20*4.04);
  double temperature=100.0;
  for(size_t i=0;i<equilibrium_time_steps/dump_inteval;i++){
    px_list.pop_front();
    py_list.pop_front();
    pz_list.pop_front();
  }
  px_temp=averagelist(px_list);
  py_temp=averagelist(py_list);
  pz_temp=averagelist(pz_list);
  for(size_t i=0;i<useful;i++){
    px_vector[i]=px_list.front()-px_temp;
    py_vector[i]=py_list.front()-py_temp;
    pz_vector[i]=pz_list.front()-pz_temp;
    px_list.pop_front();
    py_list.pop_front();
    pz_list.pop_front();
  }
  double* px_auto_corre=polarcorrelation_zero_padding(px_vector,useful);
  fftw_complex* out;
  fftw_plan p;
  out=(fftw_complex* )fftw_malloc(sizeof(fftw_complex)*((useful-1)/2+1));
  p=fftw_plan_dft_r2c_1d((useful-1),px_auto_corre,out,FFTW_ESTIMATE);
  fftw_execute(p);
  std::fstream fs_px_frequency;
  fs_px_frequency.open("px_frequency.txt",std::fstream::out);
  double first_term=dielectric_first(polarvar(px_vector,useful),volume,temperature);
  std::cout<<first_term<<std::endl;
  for(size_t i=0;i<useful/2+1;i++){
    fs_px_frequency<<out[i][0]<<" "<<out[i][1]<<std::endl;
    std::cout<<std::complex<double>(first_term,0.0)+dielectric_second(std::complex<double>(out[i][0],out[i][1]),volume,temperature,(i+0.0)/useful)<<std::endl;
  }
  fs_px_frequency.close();
}
