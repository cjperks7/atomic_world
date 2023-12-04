

import aurora
import os


filep = os.path.join(
    '/home/cjperks',
    'tofu_sparc/atomic_data/ADAS_PEC_files',
    'fsciortino_share/atomdb',
    'pec#ar17.dat'
    )

Te_eV = 10e3
ne_cm3 = 1.0

aurora.adf15_line_identification(
    filep, 
    Te_eV=Te_eV, 
    ne_cm3=ne_cm3, 
    )

trs = aurora.read_adf15(filep)

ind = np.where(trs['lambda [A]'] == 3.731117)

trs['type'][1]  # ioniz
trs['type'][15] # excit
trs['type'][61] # recomb

10**trs['log10 PEC fun'][1](np.log10(1e14), np.log10(10e3))
10**trs['log10 PEC fun'][15](np.log10(1e14), np.log10(10e3))
10**trs['log10 PEC fun'][61](np.log10(1e14), np.log10(10e3))