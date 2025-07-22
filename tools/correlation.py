import numpy as np
from mpi4py import MPI
import time
def OneDto3D(ind,cell):
  z=np.floor(ind/cell/cell);
  y=np.floor((ind-z*cell*cell)/cell);
  x=np.floor(ind-z*cell*cell-y*cell);
  return [int(x),int(y),int(z)]
def ThreeDto1D(indlist,cell):
  return indlist[0]+indlist[1]*cell+indlist[2]*cell*cell;
def searchlist(ind,cell,cutoffmax,celllength):
  incremax=int(np.ceil(cutoffmax/celllength));
  centerindex=OneDto3D(ind,cell);
  neilist=[];
  for i in range(-incremax,incremax):
    for j in range(-incremax,incremax):
      for k in range(-incremax,incremax):
        if i==0 and j==0 and k==0:
          continue;
        if i**2+j**2+k**2 < incremax**2:
          ind=ThreeDto1D([(centerindex[0]+i)%cell,(centerindex[1]+j)%cell,(centerindex[2]+k)%cell],cell);
          neilist.append(int(ind));
  return neilist
def correlation(index1,index2,cell,polartimeprofile,timeframe):
  [ind1x,ind1y,ind1z]=OneDto3D(index1,cell);
  [ind2x,ind2y,ind2z]=OneDto3D(index2,cell);
  corresum=0.0
  for i in range(timeframe):
    vect1=polartimeprofile[i][ind1z][ind1y][ind1x][:];
    vect2=polartimeprofile[i][ind2z][ind2y][ind2x][:];
    corr=np.dot(vect1,vect2)/np.linalg.norm(vect1)/np.linalg.norm(vect2);
    corresum=corr+corresum;
  return corresum/timeframe;
comm=MPI.COMM_WORLD;
rank=comm.Get_rank();
size=comm.Get_size();
data=np.fromfile("local_polar.bin",dtype=np.float64);
cell=20;
timeframe=len(data)/(3*cell**3);
timeframe=int(timeframe);
polartimeprofile=data.reshape((timeframe,cell,cell,cell,3));
celllength=4;
deltat=0.001*200;
cutoffmax=30.0;
print("@@@I am Entering the Loop@@@")
f=open("./CorreFileOnNode{0:d}.dat".format(rank),'w');
t1=time.time();
for ind in range(rank,cell**3,size):
  [k,j,i]=OneDto3D(ind,cell);
  neilist=searchlist(ind,cell,cutoffmax,celllength);
  for m in range(len(neilist)):
    corre=correlation(ind,neilist[m],cell,polartimeprofile,250);
    [neix,neiy,neiz]=OneDto3D(neilist[m],cell);
    distance=celllength*np.sqrt((neiz-i-np.round((neiz-i)/cell)*cell)**2+(neiy-j-np.round((neiy-j)/cell)*cell)**2+((neix-k-np.round((neix-k)/cell)*cell))**2);
    f.write("{0:10.7f} {1:10.7f}\n".format(distance,corre));
  print(k)
f.close();
t2=time.time();
if rank==0:
  print("Time Cost is := {0:10.7f}".format(t2-t1));
