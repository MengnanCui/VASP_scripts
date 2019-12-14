"""
Created on Sat Apr  7 19:09:15 2018

@author: mncui
version = 2.0
data = 2018.08.29
"""
import os
print(os.getcwd())
inputfile = input('input filename:  ') 
data = [line.strip(' ').split() for line in open(str(inputfile))]  # write ini.txt to data
posnew = open('POSCAR','w') 
print('Attention: now you need to input two coordinations which are dmax and dmin,respectively. Atoms that bettween them(dmin<=Atoms<=dmax) can move during relaxing!!!')
dmax = input('input dmax:  ')     #input atom height that you do not want to change
dmin = input('input dmin:  ') 
     
a = 7
dmaxpos = len(data)     # len of txt file

d = int(dmaxpos)     

for num in range(0, d): 
    if num < a:
        posnew.write('   '.join(data[num]) + '\n')
    elif num == a:
        posnew.write('Selective\nDirect\n')
    elif num > a:
        if float(data[num][2]) < float(dmin):       
            data[int(num)][2] = data[int(num)][2] + str(' F F F')
            posnew.writelines([str(data[num][0]),'   ',str(data[num][1]),'   ',str(data[num][2]),'   '+'\n'])
        elif float(data[num][2]) > float(dmax):
            data[int(num)][2] = data[int(num)][2] + str(' F F F')
            posnew.writelines([str(data[num][0]),'   ',str(data[num][1]),'   ',str(data[num][2]),'   '+'\n'])
        else:
            data[int(num)][2] = data[int(num)][2] + str(' T T T')
            posnew.writelines([str(data[num][0]),'   ',str(data[num][1]),'   ',str(data[num][2]),'   '+'\n'])
print('Well done')
#posnew.close()