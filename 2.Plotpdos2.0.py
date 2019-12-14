import os
import matplotlib.pyplot as plt
import numpy as np
import random

def s_orbit(filename):
    x = np.array([float(l.split()[0]) for l in open (str(filename))])
    Sup = [float(l.split()[1]) for l in open(str(filename))]
    Sdown = [float(l.split()[2]) for l in open(str(filename))]
    plt.plot(x, Sup,linewidth = 1, color = 'g',label = 'Pup')
    plt.plot(x, Sdown,linewidth = 1, color = 'g',label = 'Pdown')

    plt.xlabel("Energy (eV)", fontsize = 14)
    plt.ylabel("DOS", fontsize = 14)
    plt.xlim([-10, 5])
    plt.ylim([-0.75,0.75])
    plt.vlines(0,-0.75,0.75,colors = 'black',linestyle="--")
    plt.savefig('f'+str(filename) + '.png')
    plt.close()
    print("Well done!")
def sp_orbit(filename):
    x = [float(l.split()[0]) for l in open (str(filename))]
    Sup = [float(l.split()[1]) for l in open(str(filename))]
    Sdown = [float(l.split()[2]) for l in open(str(filename))]
    Pup = np.array([float(l.split()[3]) for l in open (str(filename))]) + np.array([float(l.split()[5]) for l in open (str(filename))])
 + np.array([float(l.split()[7]) for l in open (str(filename))])
    Pdown = np.array([float(l.split()[4]) for l in open (str(filename))]) + np.array([float(l.split()[6]) for l in open (str(filename))
]) + np.array([float(l.split()[8]) for l in open (str(filename))])
    s = plt.plot(x, Sup,linewidth = 1, color = 'g',label = 'Pup')
    plt.plot(x, Sdown,linewidth = 1, color = 'g',label = 'Pdown')
    p = plt.plot(x, Pup,linewidth = 1, color = 'b',label = 'Pup')
    plt.plot(x, Pdown,linewidth = 1, color = 'b',label = 'Pdown')

    plt.xlabel("Energy (eV)", fontsize = 14)
    plt.ylabel("DOS", fontsize = 14)

    plt.xlim([-10, 5])
    plt.ylim([-0.75,0.75])
    plt.vlines(0,-0.75,0.75,colors = 'black',linestyle="--")
#    plt.legend([s,p],('S','P'),'best', numpoints=1)
    plt.savefig('f'+str(filename) + '.png')
    plt.close()
    print("Well done!")
def spd_orbit(filename):
    x = [float(l.split()[0]) for l in open (str(filename))]
    Sup = [float(l.split()[1]) for l in open(str(filename))]
    Sdown = [float(l.split()[2]) for l in open(str(filename))]

    Pup = np.array([float(l.split()[3]) for l in open (str(filename))]) + np.array([float(l.split()[5]) for l in open (str(filename))])
 + np.array([float(l.split()[7]) for l in open (str(filename))])
    Pdown = np.array([float(l.split()[4]) for l in open (str(filename))]) + np.array([float(l.split()[6]) for l in open (str(filename))
]) + np.array([float(l.split()[8]) for l in open (str(filename))])

    Dup = np.array([float(l.split()[9]) for l in open (str(filename))]) + np.array([float(l.split()[11]) for l in open (str(filename))]
) + np.array([float(l.split()[13]) for l in open (str(filename))]) + np.array([float(l.split()[15]) for l in open (str(filename))]) + n
p.array([float(l.split()[17]) for l in open (str(filename))])
    Ddown = np.array([float(l.split()[10]) for l in open (str(filename))]) + np.array([float(l.split()[12]) for l in open (str(filename
))]) + np.array([float(l.split()[14]) for l in open (str(filename))]) + np.array([float(l.split()[16]) for l in open (str(filename))]) 
+ np.array([float(l.split()[18]) for l in open (str(filename))])

    plt.plot(x, Sup,linewidth = 1, color = 'g',label = 'Pup')
    plt.plot(x, Sdown,linewidth = 1, color = 'g',label = 'Pdown')
    plt.plot(x, Pup,linewidth = 1, color = 'b',label = 'Pup')
    plt.plot(x, Pdown,linewidth = 1, color = 'b',label = 'Pdown')
    plt.plot(x, Dup,linewidth = 1, color = 'r',label = 'Pup')
    plt.plot(x, Ddown,linewidth = 1, color = 'r',label = 'Pdown')

    plt.xlim([-10, 5])
    plt.ylim([-1.5,1.5])
    plt.xlabel("Energy (eV)", fontsize = 14)
    plt.ylabel("DOS", fontsize = 14)
    plt.vlines(0,-1.5,1.5,colors = 'black',linestyle="--")
    plt.savefig('f'+str(filename) + '.png')
    plt.close()
    print("Well done!")

#-------------------------------------------------------------------------------------------
print('#---------------------------------------------------------')
begin = input('begin from which atom: ')
end = input('how many atom do you want to plot: ')
orbit = input('plot which orbit(s,sp or spd): ')
filename = []
for i in range(int(end)):
    num = int(begin) + i
    dos = 'DOS'+str(num)
    filename.append(dos)
    print('ploting'+dos)
    if orbit == 's':
        s_orbit(filename[i])
    elif orbit == 'sp':
        sp_orbit(filename[i])
    elif orbit == 'spd':
        spd_orbit(filename[i])
    else:
        print('wow,you input wrong orbit!(s,sp or spb)')
