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

