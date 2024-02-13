'''

This script is meant to test using ColRadPy to read FAc atomic
data and to compare against Adam Foster's Ar PEC data

cjperks
Jan 8, 2023

'''

# Modules
import sys, pdb
path2_colradpy = '/home/cjperks/ColRadPy'
sys.path.append(path2_colradpy)
from colradpy import colradpy

%reload_ext autoreload
%autoreload 2

# Common FAC atomic data files name
file = '/home/cjperks/2008_SPARCxray/FAC/ATOM/Ar/Ar01'

# Simulation grids
temp_arr = np.logspace(np.log10(8.62), np.log10(8.5e4), 41) # [eV]
dens_arr = np.logspace(np.log10(1e10), np.log10(1e20), 21) # [cm^-3]
meta_arr = np.array([0])

# Load atomic data
Ar = colradpy(
    file,
    meta_arr,
    temp_arr,
    dens_arr,
    use_recombination=False,
    use_recombination_three_body = False,
    use_ionization=False,
    suppliment_with_ecip=False,
    # New material
    atomic_data_type = 'FAC',
    ele = 'Ar',
    nele = 1,
    Zele = 18,
    EEDF = 'Maxwellian',
    atomic_physics = 'incl_all',
    )

# Runs collisional-radiative model
Ar.solve_cr()

# Write PEC data file
Ar.write_pecs_adf15(
    fil_name='output/kr'+str(cs)+'_test.dat'
    )