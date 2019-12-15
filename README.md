# VASP_scripts
Scripts for easily using VASP
Python >= 3.0 is required. Just copy the 'name.py' into your work directory and run! No external dependencies.
## 1. Fixing.py
- This scripts was used to fix atoms' coordinates in POSCAR

```
$ python fixing.py
> input filename:   # input your file name
> Attention: now you need to input two coordinations which are dmax and dmin,respectively. Atoms that bettween them(dmin<=Atoms<=dmax) can move during relaxing!!!
> input dmax:		# Anstrom(eg:0.5)
> input dmin:		# Anstrom(eg:0.1)
> Well done			# output 'POSCAR' file
```

## 2. Plotpdos2.0.py 
- Batch plot PDOS of each atom.
- First step:  
Calculated PDOS by VASP to get 'DOSCAR'
- Second step:  
`$ split_dos`  
you can get 'DOS0,DOS1,DOS2...'
- Third step:  
```
$ python 2.Plotpdos.py
> #---------------------------------------------------------
> begin from which atom:                # eg:1
> how many atom do you want to plot:    # eg:10
> plot which orbit(s,sp or spd):        # eg:sp
> plotingDOS1
> Well done!
> plotingDOS2
> Well done!
...
> plotingDOS10
> Well done!
```
- Finally step:  
you can get "fDOS1.png,fDOS2.png,...fDOS10.png"  
`$ eog fDOS10.png`  
![PDOS image](https://github.com/mnTusi/VASP_scripts/blob/master/image.png)
## 3. Generate slab structure and high-throughput calculations
1. Preparation  
    (1) Input folder  
    - INCAR_relax   
    - vasp.new <br>

    (2) Origin folder   
    - structure1.vasp
    - structure2.vasp  
1. Begin calculation  
`$ bsub<main.lsf`
After above finishd
`./calculate.py`
3. Introduction
    - main.lsf  
This file was used to submit `surslab.1.py` into supercomputer
4. surslab.1.py
This scripts were used to generate slabs and submit every slab structure into supercomputer to relax it.
    - Calculate bulk energy
    - Generate slabs  
`slabs = generate_all_slabs(struct,2,15,15)`  
You can edit this line to control how many slbs you want to calculate.
    - Calculated surface energy of each slabs
5. calculate.py
After finished all calculation, you can run this scripts for generating file, which include the miller index of each slab and their surface enrgy.


