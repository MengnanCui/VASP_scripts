#!/gpfs/home/mncui/soft/anaconda3/bin/python3
# Import the neccesary tools to generate surfaces
from pymatgen.core.surface import SlabGenerator, generate_all_slabs, Structure, Lattice, get_d
from pymatgen.io.vasp.inputs import Incar, Kpoints, Kpoints_supported_modes, Poscar, Potcar, VaspInput
from pymatgen.io.vasp.outputs import Dynmat, Outcar, Oszicar
# Import the neccesary tools for making a Wulff shape

import os
import shutil
import re
import time
import subprocess
import pprint
f = [line.strip('\n').split() for line in open('dir','r')]

energy = open('energy','w')
for i in range(len(f)):
	path = str(f[i][0])
	energy.write(str(path)+'\t')
	os.chdir(path)
	try:
		os.system('chmod 700 cal_surface_energy.py')
		os.system('./cal_surface_energy.py')
		eng = subprocess.check_output('grep E_surf surface_energy |tail -n 1',shell=True)
		energy.write('{} \n'.format(eng))
	except IOError:
		energy.write('{} \n'.format(int(999)))
#       os.system('grep E_surf surface_energy')
	os.chdir('../../../')

