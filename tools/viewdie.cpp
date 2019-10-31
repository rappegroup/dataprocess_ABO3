#include <stdlib.h>
#include <stdio.h>
#include <iostream>
int main(){
  FILE* fp;
  fp=fopen("local_die.bin","rb");
  double localdie[3]={0.0,0.0,0.0};
  int cell=20;
  for(size_t i=0;i<cell*cell*cell;i++){
  fread(localdie,sizeof(double),3,fp);
  std::cout<<localdie[0]<<" "<<localdie[1]<<" "<<localdie[2]<<std::endl;
  }
}
