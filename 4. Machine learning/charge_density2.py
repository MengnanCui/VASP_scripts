#!/gpfs/home/mncui/soft/anaconda3/bin/python3
# coding: utf-8

# In[1]:


import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
from collections import Iterable
import math
import time


# In[71]:


def cutoff(distance):
    fc = 1/(1+np.exp(10*(distance-9)))
    return fc

def Scale(sigma,dist):
    G = np.exp((-dist**2)/(2*sigma**2))
#    plt.plot(dist,G)
#    plt.show()
    S = (G*cutoff(dist)).sum()*(1/((math.sqrt(2*math.pi)*sigma)**3))
    return S

def Vector(sigma,dist,dxyz):
    G = (np.exp((-dist**2)/(2*sigma**2)))*(dxyz/(2*sigma**2))
    V = (G*cutoff(dist)).sum()*(1/((math.sqrt(2*math.pi)*sigma)**3))
    return V


def Tensor(sigma,dist,dx,dy):
    G = (np.exp((-dist**2)/(2*sigma**2)))*(dx*dy/(4*sigma**4))
    V = (G*cutoff(dist)).sum()*(1/((math.sqrt(2*math.pi)*sigma)**3))
    return V
def inputdata(distance, sigma):
#    os.system('neighbors.pl POSCAR 127')
#    print(os.listdir())
    data = distance
    dist = data[1:,3]
#    n = data[:,0]
    symple = []
    dx = data[1:,0]-data[0,0]
    dy = data[1:,1]-data[0,1]
    dz = data[1:,2]-data[0,2]
#    f = open('input','w')
    for s in sigma:
        V= math.sqrt(Vector(s,dist,dx)**2 
                 + Vector(s,dist,dy)**2 
                 +  Vector(s,dist,dz)**2)
        Ta = Tensor(s,dist,dx,dx)**2 
        + Tensor(s,dist,dy,dy)**2 
        + Tensor(s,dist,dz,dz)**2
        Tb = Tensor(s,dist,dx,dx)*Tensor(s,dist,dy,dy)
        + Tensor(s,dist,dy,dy)*Tensor(s,dist,dz,dz)
        + Tensor(s,dist,dx,dx)*Tensor(s,dist,dz,dz)
        - Tensor(s,dist,dx,dy)*Tensor(s,dist,dx,dy)
        - Tensor(s,dist,dy,dz)*Tensor(s,dist,dy,dz)
        - Tensor(s,dist,dx,dz)*Tensor(s,dist,dx,dz)
        Mt = np.array([[Tensor(s,dist,dx,dx),Tensor(s,dist,dx,dy),Tensor(s,dist,dx,dz)],
             [Tensor(s,dist,dy,dx),Tensor(s,dist,dy,dy),Tensor(s,dist,dy,dz)],
             [Tensor(s,dist,dz,dx),Tensor(s,dist,dz,dy),Tensor(s,dist,dz,dz)]], 'float64')  
        Tc = np.linalg.det(Mt)
        with open('symple.dat','a') as f:
            f.write("{:10f}\t{:10f}\t{:10f}\t{:10f}\t{:10f}\t".format(Scale(s,dist),V,Ta,Tb,Tc))
            f.close()
    with open('symple.dat','a') as f:
        f.write('\n')
        f.close()


# In[5]:


def fft_grid():
    coor = open('coor.dat','w')
    for i in range(NGZF):
        Nz = i+1
        for h in range(NGYF):
            Ny = h+1
            for k in range(NGXF):
                Nx = k+1
                Rx = np.linalg.norm((Nx/NGXF)*a)
                Ry = np.linalg.norm((Ny/NGYF)*b)
                Rz = np.linalg.norm((Nz/NGZF)*c)
                coor.write("{:10f} \t{:10f} \t{:10f} \n".format(Rx,Ry,Rz))
    coor.close()


# In[69]:


def distance(sigma):
    grid = np.loadtxt('coor.dat')
    for gridpoint in grid[:,]:
        distance = np.r_[gridpoint,[0]].reshape((1,4))
        for i in range(int(natom)):
            lines = i+7
            lattic = np.array(posplit[lines],float).reshape((1,3))   # atom coordinate            
            dist = np.linalg.norm(lattic-gridpoint)                  # dist between grid and atom i
            tot = np.c_[lattic,dist]                                 # x y z dist
            distance = np.r_[distance,tot]
        inputdata(distance,sigma)
            
            
#            distance.write("{}\t{:10f}\t{}\n".format(i+1,dist))
#        distance.close()


# In[68]:


#-----------------------get data----------------------
log = open('run.log','w')
start_time = time.time()
log.write('obtainning data-------------------'+'\n'+'Begin: '+time.asctime( time.localtime(time.time()) )+'\n')
log.flush()
os.system('rm pos.vasp')
os.system('ase convert POSCAR pos.vasp')
posplit = [line.strip(' ').split() for line in open('pos.vasp')]
chgsplit = [line.strip(' ').split() for line in open('CHGCAR')]
# latice vector 
a, b, c = np.array(posplit[2:5],float)
# atom number
natom = sum([item for item in [float(i) for i in posplit[5][:]]])


# In[4]:


# initail line of charge density, be carefully!
ini = natom + 8 + 2
NGXF, NGYF, NGZF = chgsplit[int(ini-1)]
NGXF,NGYF,NGZF=int(NGXF),int(NGYF),int(NGZF)
log.write('End:   '+time.asctime( time.localtime(time.time()) )+'\n')
log.write('costed {} seconds'.format(time.time() - start_time)+'\n')
log.flush()

# In[6]:


# writting fft_grid-point,"coor.dat"
start_time = time.time()
log.write('writting fft_gridpoint corr.dat-------------------'+'\n'+'Begin: '+time.asctime( time.localtime(time.time()) )+'\n')
log.flush()
fft_grid()
log.write('End:   '+time.asctime( time.localtime(time.time()) )+'\n')
log.write('costed {} seconds'.format(time.time() - start_time)+'\n')
log.flush()


# In[72]:


# writting "symple.dat"
log.write('writting symple.dat-------------------'+'\n'+'Begin: '+time.asctime( time.localtime(time.time()) )+'\n')
log.flush()
sigma = np.arange(1,11,1)
distance(sigma)
log.write('End:   '+time.asctime( time.localtime(time.time()) )+'\n')
log.write('costed {} seconds'.format(time.time() - start_time)+'\n')
log.flush()


# In[ ]:


# reading and writting charge density file"chg.csv"
start_time = time.time()
log.write('writting charge density file chg.csv-------------------'+'\n'+'Begin: '+time.asctime( time.localtime(time.time()) ))
log.flush()
chgcar = open("chgcar.dat",'w')
NGF = NGXF*NGYF*NGZF
ini = 10+natom
end = (ini+math.ceil(NGF/5))
chg = np.array(chgsplit[int(ini):int(end)])
for list in chg:
    for num in list:
        chgcar.write("{:10f}\n".format(float(num)/NGF))
chgcar.close()

#chgcar = pd.DataFrame(chg)
#chgcar.to_csv('chg.csv',index =None,header=None)
log.write('\n'+'End:   '+time.asctime( time.localtime(time.time()) )+'\n')
log.write('costed {} seconds'.format(time.time() - start_time)+'\n')
log.flush()

