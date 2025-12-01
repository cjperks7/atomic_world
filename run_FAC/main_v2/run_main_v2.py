'''

Script to make FAC cross-section/rate data using
MIT Engaging

cjperks
May 10, 2024

'''

# Modules
import sys, os
from datetime import date
sys.path.insert(0, '/home/cjperks')
from atomic_world.run_FAC.main_v2 import _main as FAC
sys.path.pop(0)

nele = int(sys.argv[1])
sp = sys.argv[2]
fol = sys.argv[3]
#nele = 13
#sp = 'W'
#fol = '/home/cjperks/test'


# Prints console stdout to file
sys.stdout = open(
    os.path.join(
        fol,
        '%s%02d.out'%(sp, nele)
        ),
    'w'
    )
print('Created by: cjperks')
print('Used default settings in atomic_world')
today = date.today().strftime("%b-%d-%Y")
print("Today's date:", today)
print('Inputs:')
print(nele)
print(sp)
print(fol)
print('Outputs:')

# runs FAC
FAC.rates(
    ion=sp,
    save_path = fol,
    nele = nele,
    #settings = str(nele),
    #settings = 'EBIT_W_13',
    #settings = 'XEUS_Mo_'+str(nele),
    #settings = 'Survey_W_'+str(nele),
    #settings = 'EBIT_W_'+str(nele),
    settings = 'W7X_W_'+str(nele),
    #physics = ['en', 'ce'],
    )
