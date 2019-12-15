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

def write_poscar(stru,hkl,path):
	os.makedirs(str(path)) # make directory for 'structure/miller/poscar'
	slab.make_supercell([[2,0,0],[0,2,0],[0,0,1]])
	slab.to('poscar',str(path)+'/POSCAR')

def write_kpoints(structure,kppvol,force_gamma,path):
	# Returns an automatic Kpoint object based on a structure and a kpoint density per inverse Angstrom^3 of reciprocal cell.
	# structure-Input sturcture
	# kppvol-Grid density per Angstrom^(-3) of reciprocal cell
	# force_gamma (bool)  - Force a gamma centered mesh

	#kpoints = Kpoints.automatic_density_by_vol(structure,int(kppvol),force_gamma)
	kpoints = Kpoints.gamma_automatic(kpts=(1,1,1),shift=(0,0,0))
	kpoints.write_file(str(path)+'/KPOINTS')
	log.write('.  ')

def write_incar(caltype,path):
	# caltype - calculation type: 'relax', 'static','bader'or 'pdos'
	incar = Incar.from_file('input/INCAR'+'_'+str(caltype))
	incar = incar + {'LDIPOL':'.TRUE.','IDIPOL':'3'}
	incar.write_file(str(path)+'/INCAR')
	log.write('.  ')

def write_potcar(num,element,path):
	pot = open(str(path) + '/POTCAR','w')
	for i in range(num):
		potcar = Potcar.from_file('/gpfs/home/mncui/bin/psudopotential/paw_pbe/' + str(psudo[element[i]]) + '/POTCAR')
		pot.write(str(potcar))
	pot.close()
	log.write('.\n')

def vasp_run(path):
	shutil.copy(str(cwd)+'/input/vasp.new',str(path))
	out = subprocess.check_output('bjobs',shell=True)
	put = []
	put = re.findall(r"PEND",str(out))
	while put != []:# when jobs pending
		if len(put) < 2: # if num of pending less than 3,bsub continuum
			os.chdir(path)
			bsub = subprocess.check_output('bsub < vasp.new',shell=True)
			log.write('\n'+str(bsub)+'\n')
			log.write(time.asctime( time.localtime(time.time()) ) +'\n')
			log.flush()
			os.chdir('../../../')
			return
		else:
			log.write('Pending: >%s> ' %len(put))
			time.sleep(5)
			out = subprocess.check_output('bjobs',shell=True)
			put = []
			put = re.findall(r"PEND",str(out))
			log.write(' >%s> ' %len(put))
			log.flush()
	os.chdir(path)
	bsub = subprocess.check_output('bsub < vasp.new',shell=True)
	print(bsub)
	log.write('\n'+str(bsub)+'\n')
	log.write(time.asctime( time.localtime(time.time()) ) +'\n')
	log.flush()
	os.chdir('../../../')
psudo = {'Ac': 'Ac', 'Ag': 'Ag', 'Al': 'Al', 'Ar': 'Ar', 'As': 'As', 'Au': 'Au', 'B': 'B', 'Ba': 'Ba_sv', 'Be': 'Be_sv', 'Bi': 'Bi', 'Br': 'Br', 'C': 'C', 'Ca': 'Ca_sv', 'Cd': 'Cd', 'Ce': 'Ce', 'Cl': 'Cl', 'Co': 'Co', 'Cr': 'Cr_pv', 'Cs': 'Cs_sv', 'Cu': 'Cu_pv', 'Dy': 'Dy_3', 'Er': 'Er_3', 'Eu': 'Eu', 'F': 'F', 'Fe': 'Fe_pv', 'Ga': 'Ga_d', 'Gd': 'Gd', 'Ge': 'Ge_d', 'H': 'H', 'He': 'He', 'Hf': 'Hf_pv', 'Hg': 'Hg', 'Ho': 'Ho_3', 'I': 'I', 'In': 'In_d', 'Ir': 'Ir', 'K': 'K_sv', 'Kr': 'Kr', 'La': 'La', 'Li': 'Li_sv', 'Lu': 'Lu_3', 'Mg': 'Mg_pv', 'Mn': 'Mn_pv', 'Mo': 'Mo_pv', 'N': 'N', 'Na': 'Na_pv', 'Nb': 'Nb_pv', 'Nd': 'Nd_3', 'Ne': 'Ne', 'Ni': 'Ni_pv', 'Np': 'Np', 'O': 'O', 'Os': 'Os_pv', 'P': 'P', 'Pa': 'Pa', 'Pb': 'Pb_d', 'Pd': 'Pd', 'Pm': 'Pm_3', 'Pr': 'Pr_3', 'Pt': 'Pt', 'Pu': 'Pu', 'Rb': 'Rb_sv', 'Re': 'Re_pv', 'Rh': 'Rh_pv', 'Ru': 'Ru_pv', 'S': 'S', 'Sb': 'Sb', 'Sc': 'Sc_sv', 'Se': 'Se', 'Si': 'Si', 'Sm': 'Sm_3', 'Sn': 'Sn_d', 'Sr': 'Sr_sv', 'Ta': 'Ta_pv', 'Tb': 'Tb_3', 'Tc': 'Tc_pv', 'Te': 'Te', 'Th': 'Th', 'Ti': 'Ti_pv', 'Tl': 'Tl_d', 'Tm': 'Tm_3', 'U': 'U', 'V': 'V_pv', 'W': 'W_pv', 'Xe': 'Xe', 'Y': 'Y_sv', 'Yb': 'Yb_2', 'Zn': 'Zn', 'Zr': 'Zr_sv'}

def addsd(structure, midcoor, path):
	# structure - POSCAR
	# midcoor - middle of coordinates for fix
	lines = open(structure,'r').readlines()
	with open(structure,'w') as f:
		if midcoor < 0.5: # if atom all  on down side of slabs
			for num in range(len(lines)):
				if num > 7:
					splits = lines[num].split()
					if float(splits[2])< float(midcoor): # if atom low that midcoor, fixing
						f.write(splits[0]+' '+splits[1]+' '+splits[2]+' F F F ')
					else: # if atom above midcoor , allow movement
						f.write(splits[0]+' '+splits[1]+' '+splits[2]+' T T T ')
					if len(splits)>3: #If the atom had an identifer (ex: '.5 .5 .5 Mg'), add the identifier
						f.write(splits[3]+' \n')
					else:
						f.write(' \n')
				elif num < 7: # Write all lines prior to line 7 as is
					f.write(lines[num])
				elif num==7:
					f.write('Selective \n') # Add 'S' tag whil keeping
					f.write(lines[num]) # Add previous identifier
		if midcoor > 0.5: # if atom all on top side of slabs
			for num in range(len(lines)):
                                if num > 7:
                                        splits = lines[num].split()
                                        if float(splits[2]) > float(midcoor): # if atom low that midcoor, fixing
                                                f.write(splits[0]+' '+splits[1]+' '+splits[2]+' F F F ')
                                        else: # if atom above midcoor , allow movement
                                                f.write(splits[0]+' '+splits[1]+' '+splits[2]+' T T T ')
                                        if len(splits)>3: #If the atom had an identifer (ex: '.5 .5 .5 Mg'), add the identifier
                                                f.write(splits[3]+' \n')
                                        else:
                                                f.write(' \n')
                                elif num < 7: # Write all lines prior to line 7 as is
                                        f.write(lines[num])
                                elif num==7:
                                        f.write('Selective \n') # Add 'S' tag whil keeping
                                        f.write(lines[num]) # Add previous identifier

def surf_energy(structure, S_area, path):
	os.chdir(path)
	with open('cal_surface_energy.py','a') as f:
		f.write('#!/gpfs/home/mncui/soft/anaconda3/bin/python3.6 \n')
		f.write('import os \nimport shutil \nimport re\nimport time\nimport subprocess\n')
		f.write('from pymatgen.core.surface import SlabGenerator, generate_all_slabs, Structure, Lattice, get_d \nfrom pymatgen.io.vasp.inputs import Incar, Kpoints, Kpoints_supported_modes, Poscar, Potcar, VaspInput\nfrom pymatgen.io.vasp.outputs import Dynmat, Outcar, Oszicar\n')
		f.write('# writen by mncui 2019.02.06 for calculation Surface Energy\n')
		f.write("structure = '%s'\n" %structure)
		f.write("oszicar = Oszicar('../../bulk/'+ str(structure)+'/OSZICAR')\n")
		f.write("E_bulk = oszicar.final_energy\n")
		f.write("natom = %f \n" %natom)
		f.write("S_area = %f \n" %S_area)
		f.write("oszicars = Oszicar('OSZICAR')\n")
		f.write("E_slab = oszicars.final_energy\n")
		f.write("pos = [line.strip(' ').split() for line in open('POSCAR')]\n")
		f.write("natoms = sum([item for item in [float(i) for i in pos[6][:]]])\n")
		f.write("with open('surface_energy','w') as f:\n")
		f.write("	f.write('E_bulk  = %.2f \\nnatom  = %s\\n' %(E_bulk, natom))\n")
		f.write("	f.write('E_slab = %.2f  \\nnatoms = %s\\n' %(E_slab, natoms))\n")
		f.write("	f.write('S_area = %.2f\\n' %S_area)\n")
		f.write("	f.write('E_surf = (E_slab - (natoms/natom) * E_bulk)/(S_area * 2)\\n')\n")
		f.write("	f.write('E_surf = %.3f' %((float(E_slab) - (float(natoms)/float(natom)) * float(E_bulk))/(S_area * 2)))")
	os.chdir('../../../')
#----------------------------main.py_begin------------------------------


cwd = os.getcwd()
log = open('log','w')
dire = open('dir','w')
log.write(cwd + '\n')
poscar = os.listdir('origin/')
for stru in poscar:
	log.write('----------------------'+str(stru)+'-------------------------\n')
	origin = 'origin/' + str(stru)  # stru are name of POSCAR file from origin directory 
	struct = Structure.from_file(origin) # read POSCAR 
	i = 0
	j = 0
	hkl = []
	posplit = [line.strip(' ').split() for line in open(origin)]
	element_line = posplit[5][:]
	num = len(element_line)      # atomic type number of bulk
	natom = sum([item for item in [float(i) for i in posplit[6][:]]])   # atomic number of bulk
	#----------calculate bulk energy------------
	os.makedirs('slab/bulk/' + str(stru))
	struct.to('poscar', 'slab/bulk/'+ str(stru) + '/POSCAR')
	shutil.copy(str(cwd)+'/input/INCAR_relax', 'slab/bulk/' + str(stru)+'/INCAR')
	write_potcar(num, element_line, 'slab/bulk/' + str(stru))
	write_kpoints(struct, 500, True, 'slab/bulk/' + str(stru))
	vasp_run('slab/bulk/' + str(stru))

	#-------------Slab generator- --------------
	slabs = generate_all_slabs(struct,2,10,10)
	log.write("%s unique slab structures have been found for a max Miller index of 2" %len(slabs)+ '\n')
	for slab in slabs:
		i = i+1
		if str(slab.miller_index) != str(hkl):
			j += 1
			hkl = slab.miller_index
			log.write('%s'%j + str(hkl) + '\n')
			center = slab.center_of_mass
			surf_area = slab.surface_area  # Surface area of slab
			log.write('	surface_area is : ' + str(surf_area)+ '\n')
			n = 0
			path = 'slab/' + str(stru)+'/'+str(hkl[0])+str(hkl[1])+str(hkl[2])+'-'+str(n)
			dire.write(path+'\n')  # write dir file for calculate.py scripts
			write_poscar(stru, hkl, path)
			
			structure = Structure.from_file(str(path)+'/POSCAR')
			write_kpoints(structure, 500, True, path)
			write_incar('relax',path)
			# write potcar
			p = [line.strip(' ').split() for line in open(str(path)+'/POSCAR')]
			pelements = p[5][:]
			pn = len(pelements)
			write_potcar(pn, pelements, path)
				
			addsd(str(path)+'/POSCAR',center[2],path)
			surf_energy(stru,surf_area,path)
			vasp_run(path)
		else:
			j += 1
			log.write('euqal in miller_index but not in determination'+str(hkl)+ '\n')
			center = slab.center_of_mass
			surf_area = slab.surface_area
			log.write('surface_area is : ' + str(surf_area)+ '\n')
			n = n+1
			path = 'slab/' + str(stru)+'/'+str(hkl[0])+str(hkl[1])+str(hkl[2])+'-'+str(n)
			dire.write(path+'\n')
			write_poscar(stru, hkl, path)
			write_kpoints(structure, 500, True, path)
			structure = Structure.from_file(str(path)+'/POSCAR')
			write_incar('relax',path)
			# write potcar
			p = [line.strip(' ').split() for line in open(str(path)+'/POSCAR')]
			pelements = p[5][:]
			pn = len(pelements)
			write_potcar(pn, pelements, path)

			addsd(str(path)+'/POSCAR',center[2],path)
			surf_energy(stru,surf_area,path)
			vasp_run(path)
