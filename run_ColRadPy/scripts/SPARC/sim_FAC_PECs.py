'''

Script to use FAC cross-section data with ColRadPy collisional-radiative
modeling to do generate PECs

cjperks
Feb 20, 2024


'''

# Modules
import sys, os
from colradpy import colradpy
import numpy as np

# Ion
sp = 'Ar'
nele = 2
Znuc = 18

# Common FAC atomic data files name
file = os.path.join(
    '/home/cjperks/work/2008_SPARCxray',
    #'FAC/ATOM',
    #sp,
    'FAC/ATOM_editted_v1',
    sp+'_v12',
    sp+'%02d'%(nele)
    #'ColRadPy/ADF04/Xe',
    #'ssh41_cl_ca#xe43.dat'
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
    #atomic_physics = ['en', 'tr', 'ce', 'rr', 'ci']
    )

# Runs collisional-radiative model
crm.solve_cr()

# Write PEC data file
crm.write_pecs_adf15(
    fil_name=os.path.join(
        '/home/cjperks',
        '2008_SPARCxray/ColRadPy/FAC',
        'test_FAC_'+sp+'%02d.dat'%(nele)
        #'2008_SPARCxray/ColRadPy/ADF04/Xe',
        #'test_xe43.dat'
        ),
    wave_lims = [0.5, 10]
    )
'''
############### Demo PECs

import aurora
fil_AF = os.path.join(
    '/home/cjperks',
    'tofu_sparc/atomic_data',
    'ADAS_PEC_files/fsciortino_share/atomdb',
    'pec#'+sp.lower()+'%02d.dat'%(Znuc-nele)
    )


# Reads ADF15
trs_AF = aurora.read_adf15(fil_AF)

# Line of interest
#lamb_AF = 3.736532
lamb_AF = 3.731117
inds = np.where(trs_AF['lambda [A]'] == lamb_AF)

ind_i = inds[0][0]
ind_e = inds[0][1]
ind_r = inds[0][2]

# AF PEC interp func
PEC_Ly_i = trs_AF['log10 PEC fun'][ind_i]
PEC_Ly_e = trs_AF['log10 PEC fun'][ind_e]
PEC_Ly_r = trs_AF['log10 PEC fun'][ind_r]

# Evalues
val_Ly_i = 10**PEC_Ly_i.ev(np.log10(1e14), np.log10(temp_arr))
val_Ly_e = 10**PEC_Ly_e.ev(np.log10(1e14), np.log10(temp_arr))
val_Ly_r = 10**PEC_Ly_r.ev(np.log10(1e14), np.log10(temp_arr))

tr = 2

recomb = crm.data['input_file']['rates']['recomb']['recomb_excit'][tr,:]

fig,ax = plt.subplots(1,2)

ax[0].plot(
    temp_arr/1e3,
    val_Ly_i,
    label = 'ioniz, AF'
    )
ax[0].plot(
    temp_arr/1e3,
    val_Ly_e,
    label = 'excit, AF'
    )
ax[0].plot(
    temp_arr/1e3,
    val_Ly_r,
    label = 'recomb, AF'
    )

ax[0].plot(
    crm.data['user']['temp_grid']/1e3,
    #crm.data['cr_matrix']['q_ji'][1,0,:]*np.exp(-3318.191/crm.data['user']['temp_grid'])*1.06,
    crm.data['input_file']['rates']['excit']['ratec_cm3/s'][tr, :],
    label = 'excit, me'
    )

ax[0].plot(
    temp_arr/1e3,
    recomb,
    label = 'recomb, rr rate coeff'
    )


ax[0].grid('on')
ax[0].set_xlabel(r'$T_e$ [$keV$]')
ax[0].set_ylabel(r'PEC [$cm^3/s$]')
ax[0].set_title('Ly_a1')
leg = ax[0].legend()
leg.set_draggable('on')


ax[1].plot(
    temp_arr/1e3,
    val_Ly_i,
    label = 'ioniz, AF'
    )
ax[1].plot(
    temp_arr/1e3,
    val_Ly_e,
    label = 'excit, AF'
    )
ax[1].plot(
    temp_arr/1e3,
    val_Ly_r,
    label = 'recomb, AF'
    )

ax[1].plot(
    crm.data['user']['temp_grid']/1e3,
    #crm.data['cr_matrix']['q_ji'][1,0,:]*np.exp(-3318.191/crm.data['user']['temp_grid'])*1.06,
    crm.data['input_file']['rates']['excit']['ratec_cm3/s'][tr, :]*1.06,
    label = 'excit, me'
    )

ax[1].plot(
    temp_arr/1e3,
    recomb,
    label = 'recomb, rr rate coeff'
    )


ax[1].grid('on')
ax[1].set_xlabel(r'$T_e$ [$keV$]')
ax[1].set_ylabel(r'PEC [$cm^3/s$]')
ax[1].set_title('Ly_a1')

ax[1].set_yscale('log')
ax[1].set_xscale('log')
ax[1].set_ylim(1e-20,1e-10)



upr = 1
#plotting the populating levels
figp, axp = plt.subplots()


axp.plot(
    (
        crm.data['processed']['pop_lvl'][upr,:,0,10,0]/
        np.sum(crm.data['processed']['pop_lvl'][upr,:,0,10,0])
        ),
    label = crm.data['user']['temp_grid'][10]/1e3
    )

axp.plot(
    (
        crm.data['processed']['pop_lvl'][upr,:,0,10,0]/
        np.sum(crm.data['processed']['pop_lvl'][upr,:,0,10,0])
        ),
    label = crm.data['user']['temp_grid'][20]/1e3
    )

axp.plot(
    (
        crm.data['processed']['pop_lvl'][upr,:,0,10,0]/
        np.sum(crm.data['processed']['pop_lvl'][upr,:,0,10,0])
        ),
    label = crm.data['user']['temp_grid'][30]/1e3
    )

axp.plot(
    (
        crm.data['processed']['pop_lvl'][upr,:,0,-1,0]/
        np.sum(crm.data['processed']['pop_lvl'][upr,:,0,-1,0])
        ),
    label = crm.data['user']['temp_grid'][-1]/1e3
    )
axp.legend()
#axp.set_xlim(0,20)






fig2, ax2= plt.subplots()

ax2.plot(
    crm.data['user']['temp_grid']/1e3,
    crm.data['input_file']['rates']['excit']['ratec_cm3/s'][0,:],
    label = 'P0.5 to grd'
    )
ax2.plot(
    crm.data['user']['temp_grid']/1e3,
    crm.data['input_file']['rates']['excit']['ratec_cm3/s'][2,:]/2,
    label = 'P1.5 to grd'
    )

tr = 2
crm.data['input_file']['rates']['excit']['ratec_cm3/s'][tr,:]/crm.data['input_file']['rates']['a_val'][tr]


fig3, ax3 = plt.subplots()

ax3.plot(
    crm.data['user']['temp_grid'],
    crm.data['processed']['pops'][0,0,:,8],
    label = '
)


from pfac import rfac
rr = rfac.read_rr(file+'a.rr')

import scipy.constants as cnt
# Useful constants
fine_struct2 = cnt.physical_constants['fine-structure constant'][0]**2 # []
eV2Hartree = 1/cnt.physical_constants['Hartree energy in eV'][0] # [Hartree/eV]
a02 = cnt.physical_constants['Bohr radius'][0]**2 *1e4 # [cm^2]

blk = 1

Einc_FAC = rr[1][blk]['EGRID']
dE = rr[1][blk]['Delta E'][0]
ion_L = 0
w_upr = rr[1][blk]['free_2J'][0]/2
w_lwr = rr[1][blk]['bound_2J'][0]/2
limit = rr[1][blk]['parameters'][0,:]

E_gamma = np.logspace(
    np.log10(dE),
    np.log10(dE*1e3)
    )
E_e = E_gamma - dE

# Formula for bound-free oscillator strength, [1/Hartree]
xx = (E_e+ limit[3]) /limit[3]
yy = (1 +limit[2]) /(np.sqrt(xx) + limit[2])

dgf_dE = (
    E_gamma /(E_e+limit[3])
    * limit[0]
    * xx**(-3.5 -ion_L +0.5 *limit[1])
    * yy**limit[1]
    )

eps = E_e * eV2Hartree
omega = E_gamma *eV2Hartree

# Function for photo-ionization cross-section, [atomic units]
XS_PI = (
    2 *np.pi *np.sqrt(fine_struct2)
    /(1 +2*w_lwr)
    *(1 +fine_struct2 *eps)
    /(1 +0.5 *fine_struct2 *eps)
    *dgf_dE
    )

# Function for radiative recombinarion cross-section, [cm2]
XS_RR = (
    fine_struct2/2
    * (1 +2*w_lwr)/(1 +2*w_upr)
    * omega**2
    / eps
    / (1 +0.5*fine_struct2 *eps)
    * XS_PI
    ) *a02



import colradpy.read_FAC_utils as utils
XSdata = utils._read_ascii(
    fil = file,
    react = 'rr',
    )

#lbl = 'ioniz'
lbl = 'recomb'
ntr = 1
ntr2 = 1

plt.loglog(
    crm.data['input_file']['rates'][lbl]['engy_eV'][ntr,:,0],
    crm.data['input_file']['rates'][lbl]['XS_cm2'][ntr,:, 0],
    '*-'
    )

plt.loglog(
    crm.data['input_file']['rates'][lbl]['engy_eV'][ntr,:,-1],
    crm.data['input_file']['rates'][lbl]['XS_cm2'][ntr,:, -1],
    '*-'
    )

plt.loglog(
    XSdata[ntr2][4]['engy'],
    XSdata[ntr2][4]['XS'],
    '*-'
    )

plt.loglog(
    E_e,
    XS_RR,
    'r'
)


import colradpy.read_FAC_utils as utils
XSdata = utils._read_ascii(
    fil = file,
    react = 'ce',
    )

#lbl = 'ioniz'
#lbl = 'recomb'
lbl = 'excit'
ntr = 2
ntr2 = 3

plt.loglog(
    crm.data['input_file']['rates'][lbl]['engy_eV'][ntr,:,0],
    crm.data['input_file']['rates'][lbl]['XS_cm2'][ntr,:, 0],
    '*-'
    )

plt.loglog(
    crm.data['input_file']['rates'][lbl]['engy_eV'][ntr,:,-1],
    crm.data['input_file']['rates'][lbl]['XS_cm2'][ntr,:, -1],
    '*-'
    )

plt.loglog(
    XSdata[0][ntr2]['engy'],
    XSdata[0][ntr2]['XS'],
    '*-'
    )
'''