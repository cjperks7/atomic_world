'''

Test script for running collisional-radiatve modeling

cjperks

Dec. 11th, 2023

TO DO:
    1) Collisional ionization (rt6)
    2) Radiative recombination (rt4)
    3) Dielectronic capture (rt5)
    4) Higher energy levels beyond nmax=2 (cascades)
    5) Transport modeling in ionization balance equation
    6) Solve many charge states at once
    7) Couple to pfac.crm to get rates
    8) Compare to Adam Foster's results

'''

# Modules
import aurora
from atomic_world.run_FAC.main import _crm as crm
import os

# Enables automatic reloading of modules
%reload_ext autoreload
%autoreload 2

# Plasma parameters
ne_cm3 = 1e14
Te_eV = 1370
ni_cm3 = 1e10

ion = 'Ar'
ZZ = 17

# Path to FAC results
path_fac = '/home/cjperks/atomic_world/run_FAC/scripts/Ar/test_ai'

# Path to adf15
path_adas = os.path.join(
    '/home/cjperks/tofu_sparc',
    'atomic_data/ADAS_PEC_files',
    'fsciortino_share/atomdb',
    'pec#ar17.dat'
    )

###################################################
#
#           Ionization Balance Modeling
#
###################################################

# States with no transport, so ionization balance
atom_data = aurora.get_atom_data(ion,['scd','acd'])

_, fz = aurora.get_frac_abundances(
    atom_data, 
    ne_cm3, 
    Te_eV, 
    #rho=rhop, 
    plot=False
    )

# Charge-state resolved ion density
nz_cm3 = ni_cm3 * fz[0,:] # dim(cs,)

#nz_cm3 = np.ones(len(fz[0,:]))
#nz_cm3[18] = 0.0
#nz_cm3[16] = 0.0

###################################################
#
#           Popullation Balance Modeling
#
###################################################

nsim = 2

# Formats FAC data
data, lvls = crm.prep(
    ion = ion,
    # FAC controls
    path_fac = path_fac,
    nele = len(nz_cm3) -1 -ZZ,
    nsim = nsim,  # Maximum quantum number to solve for population
    mode = 'mr',
    Te_eV = Te_eV,
    )

# Solves relative population balance
pops, mat, rhs = crm.pops(
    data = data,
    lvls = lvls,    # [list], levels to solve for
    nsim = nsim,    
    # Plasma
    ne_cm3 = ne_cm3,
    nz_cm3 = fz[0,:],  # dim(cs,)
    ZZ = ZZ,
    ) # [1/cm3]

# Finds PEC terms, [cm3/s]
exc, rec, ion = crm.PECs(
    # States for transition of interest
    state = 3,
    final = 0,
    # Data
    data = data,
    pops = pops,
    nsim = nsim,   
    # Plasma
    ne_cm3 = ne_cm3,
    nz_cm3 = fz[0,:],
    ZZ = ZZ,
    )


###################################################
#
#         Compare data
#
###################################################

# Read adf15 file
trs = aurora.read_adf15(path_adas)

# Wavelength of interest
lamb0 = 3.731117 # [A]

# Transitions of interest
sel = np.where(trs['lambda [A]'] == lamb0)[0]

# PECs indices
ind_i = sel[0]
ind_e = sel[1]
ind_r = sel[2]

# PEC vaues
PEC_i = 10**trs['log10 PEC fun'][ind_i](np.log10(ne_cm3), np.log10(Te_eV))
PEC_e = 10**trs['log10 PEC fun'][ind_e](np.log10(ne_cm3), np.log10(Te_eV))
PEC_r = 10**trs['log10 PEC fun'][ind_r](np.log10(ne_cm3), np.log10(Te_eV))

