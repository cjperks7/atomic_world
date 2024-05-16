'''

Script to use FAC cross-section data with ColRadPy collisional-radiative
modeling to do generate PECs

cjperks
Feb 20, 2024


'''

# Modules
import sys, os
import numpy as np
from colradpy import colradpy

# Ion
nele = int(sys.argv[1])
sp = sys.argv[2]
fol = sys.argv[3]

if sp == 'Ar':
    Znuc = 18
    lims = [3.5, 4.5]
elif sp == 'Kr':
    Znuc = 36
    lims = [0.5, 1.5]
elif sp == 'Xe':
    Znuc = 54
    if nele <= 3:
        lims = [2.0, 2.5]
    else:
        lims = [2.6, 2.8]
elif sp == 'W':
    Znuc = 74


# Common FAC atomic data files name
file = os.path.join(
    '/nobackup1/cjperks/work',
    'FAC/best',
    sp,
    sp+'%02d'%(nele)
    )

# Simulation grids
temp_arr = np.logspace(np.log10(8.62), np.log10(8.5e4), 41) # [eV]
dens_arr = np.logspace(np.log10(1e10), np.log10(1e20), 21) # [cm^-3]
meta_arr = np.array([0])

# Load atomic data
crm = colradpy(
    file,
    meta_arr,
    temp_arr,
    dens_arr,
    use_recombination=True,
    use_recombination_three_body = False,
    use_ionization=True,
    suppliment_with_ecip=False,
    # New material
    atomic_data_type = 'FAC',       # Flag data source
    ele = sp,                       # Ion species
    nele = nele,                    # Number of electrons
    Zele = Znuc,                    # Nuclear charge
    EEDF = 'Maxwellian',            # Electron energy distriubtion function
    atomic_physics = 'incl_all',    # Atomic data files to search for
    #atomic_physics = ['en', 'tr', 'ce', 'rr', 'ai'],
    )

print('ColRadPy FAC read done')

# Runs collisional-radiative model
crm.solve_cr()

print('ColRadPy CR solve done')

# Write PEC data file
crm.write_pecs_adf15(
    fil_name=os.path.join(
        fol,
        'FAC_'+sp+'%02d.dat'%(nele)
        ),
    wave_lims = lims
    )

print('ColRadPy adf15 write done')

# Saves effective rate data
np.savez(
    os.path.join(
        fol,
        'FAC_'+sp+'%02d_GCR.npz'%(nele)
        ),
    acd = crm.data['processed']['acd'],
    scd = crm.data['processed']['scd'],
    plt = crm.data['processed']['plt'],
    prb = crm.data['processed']['prb'],
    Te_eV = temp_arr,
    ne_cm3 = dens_arr
    )

print('Save ColRadPy rates done')