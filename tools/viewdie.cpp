#include <stdlib.h>
#include <stdio.h>
#include <iostream>
int main(){
  FILE* fp;
  fp=fopen("local_polar.bin","rb");
  double localdie[3]={0.0,0.0,0.0};
  int cell=10;
  int frame=20000/200;
  double** sumlocal=new double* [cell*cell*cell];
  for(size_t i=0;i<cell*cell*cell;i++){
    sumlocal[i]=new double [3];
    for(size_t j=0;j<3;j++){
      sumlocal[i][j]=0.0;
    }
  }
  for(size_t t=0;t<frame;t++){
    for(size_t i=0;i<cell*cell*cell;i++){
      fread(localdie,sizeof(double),3,fp);
      for(size_t j=0;j<3;j++){
        sumlocal[i][j]=sumlocal[i][j]+localdie[j];
      }
    }
  }
  for(size_t i=0;i<cell*cell*cell;i++){
    for(size_t j=0;j<3;j++){
      sumlocal[i][j]=sumlocal[i][j]/frame;
    }
    std::cout<<sumlocal[i][0]<<" "<<sumlocal[i][1]<<" "<<sumlocal[i][2]<<" "<<std::endl;
  }
}
