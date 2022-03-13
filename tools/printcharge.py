import numpy as np
f=open("CHARGE.dat",'w');
cell=2;
# print Asite charge;
typelist=[1,2,3,4];
chargelist=[[2.7,2.7,2.7],[3.0,3.0,3.0],[6.7,6.7,6.7],[-3.2,-3.2,-3.2]]
lammpsdatainput=open("mixdata.BTO",'r');
datalines=lammpsdatainput.readlines();
typelist=[];
for i in range(len(datalines)):
    if datalines[i].find("Atoms")!=-1:
        for j in range(cell*cell*cell*5):
            lin=datalines[i+j+2];
            infolist=lin.split();
            typelist.append(int(infolist[2])-1);
lammpsdatainput.close();
for i in range(cell*cell*cell):
    charge=chargelist[typelist[i]];
    f.write("{0:10.7f} {1:10.7f} {2:10.7f}\n".format(charge[0],charge[1],charge[2]));
for i in range(cell*cell*cell,2*cell**3):
    charge=chargelist[typelist[i]];
    f.write("{0:10.7f} {1:10.7f} {2:10.7f}\n".format(charge[0],charge[1],charge[2]));
for i in range(2*cell**3,3*cell**3):
    charge=chargelist[typelist[i]];
    f.write("{0:10.7f} {1:10.7f} {2:10.7f}\n".format(charge[0],charge[1],charge[2]));
for i in range(3*cell**3,4*cell**3):
    charge=chargelist[typelist[i]];
    f.write("{0:10.7f} {1:10.7f} {2:10.7f}\n".format(charge[0],charge[1],charge[2]));
for i in range(4*cell**3,5*cell**3):
    charge=chargelist[typelist[i]];
    f.write("{0:10.7f} {1:10.7f} {2:10.7f}\n".format(charge[0],charge[1],charge[2]));
f.close();
lammpsdatainput.close();