#include "space.h"
#include <random>
#include "atom.h"
#include <math.h>
#include <iostream>
#include <map>
double* outprod(double* b,double* c){
       double* re=new double[3];
       re[0]=b[1]*c[2]-c[1]*b[2];
       re[1]=-1*(b[0]*c[2]-c[0]*b[2]);
       re[2]=b[0]*c[1]-c[0]*b[1];
      return re;
}
double inner(double* a,double* b){
   double temp=0.0;
   for(size_t i=0;i<3;i++){
      temp=temp+a[i]*b[i];
   }
   return temp;
}
double tetrahedron(double** nei,int tickA,int tickoxy1,int tickoxy2,int tickoxy3,double* p){
    double* temp=new double[3];
    for(size_t i=0;i<3;i++){
      temp[i]=0.0;
    }
    double* a=distance(temp,nei[tickoxy1]);
    double* b=distance(temp,nei[tickoxy2]);
    double* c=distance(temp,nei[tickoxy3]);
    double sum=0.0;
    sum=sum+a[0]*(b[1]*c[2]-c[1]*b[2]);
    sum=sum-a[1]*(b[0]*c[2]-c[0]*b[2]);
    sum=sum+a[2]*(b[0]*c[1]-c[0]*b[1]);
    return fabs(sum/6);
}
/*nei has already been consider as periodical boundary condition in the previous calculation*/
int convex(double** nei,double* p,int a,int b,int c){
    double* ab=distance(nei[a],nei[b]);
    double* ac=distance(nei[a],nei[c]);
    double* n=outprod(ab,ac);
    int first;
    double center[3]={0.0,0.0,0.0};
    for(size_t i=0;i<12;i++){
       center[0]=nei[i][0]+center[0];
       center[1]=nei[i][1]+center[1];
       center[2]=nei[i][2]+center[2];
    }
    center[0]/=12;
    center[1]/=12;
    center[2]/=12;
   double* temp1;
   double* temp2;
   for(int k=0;k<12;k++){
      if(k==a || k==b || k==c){
         continue;
      }
      temp1=distance(nei[a],center);
      temp2=distance(nei[a],nei[k]);
      if(inner(temp1,n)*inner(temp2,n)<0){
         return 0;
      }
   }
   return 1;
}
double Pentahedron(atom* A,atom* oxygen,int tickA,int cell,double* p){
/*decompose this Pentahedron into two tetrahedron use the longest oxygen-oxygen edge as cut*/
    int* oxy=neighbor_o_forA(tickA,cell);
    double* nei[12];
    for(int i=0;i<12;i++){
        /*only this time we need to use periodical boundary condition*/
        nei[i]=distance(A+tickA,oxygen+oxy[i],p);
    }
    double noise[12][3];
    for(size_t i=0;i<12;i++)
        for(size_t j=0;j<3;j++){
            noise[i][j]=((rand()+0.0001)/RAND_MAX-0.5)*0.01;
            nei[i][j]=nei[i][j]+noise[i][j];
        }
    int num=12*11*10/6;
    double surfall[num][3];
    int convextick[num];
    size_t temp=0;
    for(size_t i=0;i<12;i++)
        for(size_t j=0;j<i;j++)
            for(size_t k=0;k<j;k++)
            {
                surfall[temp][0]=i;
                surfall[temp][1]=j;
                surfall[temp][2]=k;
                convextick[temp]=convex(nei,p,i,j,k);
                temp++;
    }
    /*denoise */
    for(size_t i=0;i<12;i++)
       for(size_t j=0;j<3;j++)
          nei[i][j]=nei[i][j]-noise[i][j];
    /*adding volume all tegother with each tetrahedron*/
    double sum=0;
    for(size_t i=0;i<num;i++){
       if(convextick[i]==1){
         sum=sum+tetrahedron(nei,tickA,surfall[i][0],surfall[i][1],surfall[i][2],p);
       }
    }
    return sum;
}
