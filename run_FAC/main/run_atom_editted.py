'''

run_atom.py is a script meant to facilitate running the default
FAC settings stored in the ATOM class

cjperks
Dec 20th, 2023


'''

# Modules
import sys, os
from datetime import date

nele = int(sys.argv[1])
sp = sys.argv[2]
fol = sys.argv[3]

if nele == 1:
    from pfac import atom
else:
    from pfac import atom_editted as atom

# Prints console stdout to file
sys.stdout = open(
    os.path.join(
        fol,
        '%s%02d.out'%(sp, nele)
        ),
    'w'
    )
print('Created by: cjperks')
print('Used default settings in pfac.atom class')
today = date.today().strftime("%b-%d-%Y")
print("Today's date:", today)
print('Inputs:')
print(nele)
print(sp)
print(fol)
print('Outputs:')

# Runs ATOM class
atom.atomic_data(
    nele,
    sp, 
    dir=fol
    )
