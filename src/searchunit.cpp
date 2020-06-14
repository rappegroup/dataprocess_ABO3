#include "atom.h"
#include "polarconfig.h"
#include <algorithm>
#include <string>
#include <iostream>
#include <sstrream>
#include <algorithm>
#include <cmath>
void unitcell(std::string file,int cell){
  std::fstream fs;
  fs.open(file.c_str(),std::fstream::in);
  std::string temp;
  int tempint;
  double tempdouble;
  double tempdouble2;
  double period[3]={0.0,0.0,0.0};
  mapunit=new int* [cell*cell*cell];
  std::stringstream ss;
  for(size_t i=0;i<cell*cell*cell;i++){
    mapunit[i]=new int [6+8];//1 Fe with 6 Oxygen, 8 Bi;
  }
  asite=new double* [cell*cell*cell];
  bsite=new double* [cell*cell*cell];
  osite=new double* [cell*cell*cell*3];
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
    getline(fs,temp)
    for(size_t i=0;i<cell*cell*cell;i++){
      getline(fs,temp);
      asite[i]=new double[3];
      ss.str(temp);
      for(size_t j=0;j<3;j++){
        ss>>tempint;
      }
      ss>>tempdouble;
      for(size_t j=0;j<3;j++){
        ss>>asite[j];
      }
    }
    for(size_t i=0;i<cell*cell*cell;i++){
      getline(fs,temp);
      asite[i]=new double[3];
      ss.str(temp);
      for(size_t j=0;j<3;j++){
        ss>>tempint;
      }
      ss>>tempdouble;
      for(size_t j=0;j<3;j++){
        ss>>bsite[j];
      }
    }
    for(size_t i=0;i<3*cell*cell*cell;i++){
      getline(fs,temp);
      asite[i]=new double[3];
      ss.str(temp);
      for(size_t j=0;j<3;j++){
        ss>>tempint;
      }
      ss>>tempdouble;
      for(size_t j=0;j<3;j++){
        ss>>bsite[j];
      }
    }
  }
  }
  /*search the neighbors using the sort*/
  std::vector<std::pair<double,int> > sequence(cell*cell*cell,std::pair<double,int>(0.0,0));
  for(size_t i=0;i<cell*cell*cell;i++){
    for(size_t j=0;j<cell*cell*cell;j++){
      tempdobule=far(asite[i],bsite[i],period);
      sequence[j]=std::pair<double,int>(tempdouble,j);
    }

  }
}
