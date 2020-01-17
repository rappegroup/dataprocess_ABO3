#include <iostream>
#include <sstream>
#include <string>
#include <list>
#include <fstream>
#include <new>
#include <fftw3.h>
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
  for(size_t i=0;i<equilibrium_time_steps/dump_inteval;i++){
    px_list.pop_front();
    py_list.pop_front();
    pz_list.pop_front();
  }
  for(size_t i=0;i<useful;i++){
    px_vector[i]=px_list.front();
    py_vector[i]=py_list.front();
    pz_vector[i]=pz_list.front();
    px_list.pop_front();
    py_list.pop_front();
    pz_list.pop_front();
  }
  fftw_complex* out;
  fftw_plan p;
  out=(fftw_complex* )fftw_malloc(sizeof(fftw_complex)*(useful/2+1));
  p=fftw_plan_dft_r2c_1d(useful,px_vector,out,FFTW_ESTIMATE);
  fftw_execute(p);
  std::fstream fs_px_frequency;
  fs_px_frequency.open("px_frequency.txt",std::fstream::out);
  for(size_t i=0;i<useful/2+1;i++){
    fs_px_frequency<<out[i][0]/useful<<" "<<out[i][1]/useful<<std::endl;
  }
  fs_px_frequency.close();
  fftw_destroy_plan(p);
  p=fftw_plan_dft_r2c_1d(useful,py_vector,out,FFTW_ESTIMATE);
  fftw_execute(p);
  std::fstream fs_py_frequency;
  fs_py_frequency.open("py_frequency.txt",std::fstream::out);
  for(size_t i=0;i<useful/2+1;i++){
    fs_py_frequency<<out[i][0]/useful<<" "<<out[i][1]/useful<<std::endl;
  }
  fs_py_frequency.close();
  fftw_destroy_plan(p);
  p=fftw_plan_dft_r2c_1d(useful,pz_vector,out,FFTW_ESTIMATE);
  fftw_execute(p);
  std::fstream fs_pz_frequency;
  fs_pz_frequency;
  fs_pz_frequency.open("pz_frequency.txt",std::fstream::out);
  for(size_t i=0;i<useful/2+1;i++){
    fs_pz_frequency<<out[i][0]/useful<<" "<<out[i][1]/useful<<std::endl;
  }
  fs_pz_frequency.close();
  fftw_destroy_plan(p);
}
