'''

Scripts to illustrate how I grab cross-section data

cjperks
Mar 1st, 2024

'''

# Modules
from colradpy import colradpy
from pfac import rfac
import sys,os
import scipy.constants as cnt

# Useful constants
fine_struct2 = cnt.physical_constants['fine-structure constant'][0]**2 # []
eV2Hartree = 1/cnt.physical_constants['Hartree energy in eV'][0] # [Hartree/eV]
a02 = cnt.physical_constants['Bohr radius'][0]**2 *1e4 # [cm^2]

# Ion
sp = 'Ar'
nele = 1
Znuc = 18

# Common FAC atomic data files name
file = os.path.join(
    '/home/cjperks/2008_SPARCxray',
    'FAC/ATOM',
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
    atomic_data_type = 'FAC',
    ele = sp,
    nele = nele,
    Zele = Znuc,
    EEDF = 'Maxwellian',
    atomic_physics = 'incl_all',
    )


file6 = os.path.join(
    '/home/cjperks/2008_SPARCxray',
    'FAC/ATOM',
    sp+'_v6',
    sp+'%02d'%(nele)
    )

ce = rfac.read_ce(file6+'a.ce')

E_grid = crm.data['input_file']['rates']['excit']['engy_eV'][0,:,-1]

k02 = 2 * E_grid *eV2Hartree *(
    1+
    0.5 *fine_struct2 *E_grid *eV2Hartree
    ) # [atomic units]

omega = (
    ce[1][0]['bethe'][None,:]
    * np.log(E_grid[:,None]/ce[1][0]['TE0'])
    + ce[1][0]['born'][:,0][None,:]
    ) # dim(nE, ntr)

XS_tmp = (
    np.pi *omega
    / k02[:,None]
    / (1 +ce[1][0]['lower_2J'][None,:])
    ) * a02 # [cm2], dim(nE, ntr)

fig1, ax1 =plt.subplots()

ax1.loglog(
    E_grid,
    XS_tmp[:,0],
    label = '2S_0.5->2P_0.5, formula'
    )

ax1.loglog(
    crm.data['input_file']['rates']['excit']['engy_eV'][0,:,-1],
    crm.data['input_file']['rates']['excit']['XS_cm2'][0,:,-1],
    label = '2S_0.5->2P_0.5, ColRadPy'
    )

ax1.loglog(
    ce[1][0]['EGRID']+ce[1][0]['TE0'],
    ce[1][0]['crosssection'][0,:]*1e-20,
    '*',
    label = '2S_0.5->2P_0.5, FAC'
    )

ax1.loglog(
    E_grid,
    XS_tmp[:,1],
    label = '2S_0.5->2S_0.5, formula'
    )

ax1.loglog(
    crm.data['input_file']['rates']['excit']['engy_eV'][1,:,-1],
    crm.data['input_file']['rates']['excit']['XS_cm2'][1,:,-1],
    label = '2S_0.5->2S_0.5, ColRadPy'
    )

ax1.loglog(
    ce[1][0]['EGRID']+ce[1][0]['TE0'],
    ce[1][0]['crosssection'][1,:]*1e-20,
    '*',
    label = '2S_0.5->2S_0.5, FAC'
    )

ax1.loglog(
    E_grid,
    XS_tmp[:,2],
    label = '2S_0.5->2P_1.5, formula'
    )

ax1.loglog(
    crm.data['input_file']['rates']['excit']['engy_eV'][2,:,-1],
    crm.data['input_file']['rates']['excit']['XS_cm2'][2,:,-1],
    label = '2S_0.5->2P_1.5, ColRadPy'
    )

ax1.loglog(
    ce[1][0]['EGRID']+ce[1][0]['TE0'],
    ce[1][0]['crosssection'][2,:]*1e-20,
    '*',
    label = '2S_0.5->2P_1.5, FAC'
    )

ax1.set_title('H-like Ar, excit')
ax1.grid('on')
leg = ax1.legend()
leg.set_draggable('on')
ax1.set_xlim(1e3,1e7)
ax1.set_ylim(1e-24, 2e-21)
ax1.set_xlabel('E_incident [eV]')
ax1.set_ylabel('cross-section [cm2]')






rr = rfac.read_rr(file6+'a.rr')

Ee_grid = np.logspace(np.log10(1e0), np.log10(1e5), int(1e5))

E_gamma1 = Ee_grid[:,None] + rr[1][0]['Delta E'][None,:]
E_gamma2 = Ee_grid[:,None] + rr[1][1]['Delta E'][None,:]

xx1 = (
    (Ee_grid[:,None] + rr[1][0]['parameters'][None,:,3])
    / rr[1][0]['parameters'][None,:,3]
    )
xx2 = (
    (Ee_grid[:,None] + rr[1][1]['parameters'][None,:,3])
    / rr[1][1]['parameters'][None,:,3]
    )
yy1 = (
    (1+rr[1][0]['parameters'][None,:,2])
    /(np.sqrt(xx1) + rr[1][0]['parameters'][None,:,2])
    )
yy2 = (
    (1+rr[1][1]['parameters'][None,:,2])
    /(np.sqrt(xx2) + rr[1][1]['parameters'][None,:,2])
    )

dgf_dE1 = (
    E_gamma1 /(Ee_grid[:,None] +rr[1][0]['parameters'][None,:,3])
    * rr[1][0]['parameters'][None,:,0]
    * xx1**(-3.5 -0 +0.5 *rr[1][0]['parameters'][None,:,1])
    * yy1**rr[1][0]['parameters'][None,:,1]
    )
dgf_dE2 = (
    E_gamma2 /(Ee_grid[:,None] +rr[1][1]['parameters'][None,:,3])
    * rr[1][1]['parameters'][None,:,0]
    * xx2**(-3.5 -0 +0.5 *rr[1][1]['parameters'][None,:,1])
    * yy2**rr[1][1]['parameters'][None,:,1]
    )

eps = Ee_grid *eV2Hartree # [Hartree]
omega1 = E_gamma1 *eV2Hartree # [Hartree]
omega2 = E_gamma2 *eV2Hartree # [Hartree]

XS_PI1 = (
    2 *np.pi *np.sqrt(fine_struct2)
    /(1 +rr[1][0]['bound_2J'][None,:])
    *(1 +fine_struct2 *eps[:,None])
    /(1 +0.5 *fine_struct2 *eps[:,None])
    *dgf_dE1
    )
XS_PI2 = (
    2 *np.pi *np.sqrt(fine_struct2)
    /(1 +rr[1][1]['bound_2J'][None,:])
    *(1 +fine_struct2 *eps[:,None])
    /(1 +0.5 *fine_struct2 *eps[:,None])
    *dgf_dE2
    )

XS_tmp1 = (
    fine_struct2/2
    * (1 +rr[1][0]['bound_2J'][None,:])/(1 +rr[1][0]['free_2J'][None,:])
    * omega1**2
    / eps[:,None]
    / (1 +0.5*fine_struct2 *eps[:,None])
    * XS_PI1
    ) *a02
XS_tmp2 = (
    fine_struct2/2
    * (1 +rr[1][1]['bound_2J'][None,:])/(1 +rr[1][1]['free_2J'][None,:])
    * omega2**2
    / eps[:,None]
    / (1 +0.5*fine_struct2 *eps[:,None])
    * XS_PI2
    ) *a02


fig2, ax2 = plt.subplots()

ax2.loglog(
    Ee_grid,
    XS_tmp1[:,0],
    label = 'bare->(1s)2S_0.5, formula'
    )

ax2.loglog(
    crm.data['input_file']['rates']['recomb']['engy_eV'][0,:,-1],
    crm.data['input_file']['rates']['recomb']['XS_cm2'][0,:,-1],
    label = 'bare->(1s)2S_0.5, ColRadPy'
    )

ax2.loglog(
    rr[1][0]['EGRID'],
    rr[1][0]['RR crosssection'][0,:]*1e-20,
    '*',
    label = 'bare->(1s)2S_0.5, FAC'
    )

ax2.loglog(
    Ee_grid,
    XS_tmp2[:,0],
    label = 'bare->2P_0.5, formula'
    )

ax2.loglog(
    crm.data['input_file']['rates']['recomb']['engy_eV'][1,:,-1],
    crm.data['input_file']['rates']['recomb']['XS_cm2'][1,:,-1],
    label = 'bare->2P_0.5, ColRadPy'
    )

ax2.loglog(
    rr[1][1]['EGRID'],
    rr[1][1]['RR crosssection'][0,:]*1e-20,
    '*',
    label = 'bare->2P_0.5, FAC'
    )

ax2.loglog(
    Ee_grid,
    XS_tmp2[:,1],
    label = 'bare->(2s)2S_0.5, formula'
    )

ax2.loglog(
    crm.data['input_file']['rates']['recomb']['engy_eV'][2,:,-1],
    crm.data['input_file']['rates']['recomb']['XS_cm2'][2,:,-1],
    label = 'bare->(2s)2S_0.5, ColRadPy'
    )

ax2.loglog(
    rr[1][1]['EGRID'],
    rr[1][1]['RR crosssection'][1,:]*1e-20,
    '*',
    label = 'bare->(2s)2S_0.5, FAC'
    )

ax2.loglog(
    Ee_grid,
    XS_tmp2[:,2],
    label = 'bare->2P_1.5, formula'
    )

ax2.loglog(
    crm.data['input_file']['rates']['recomb']['engy_eV'][3,:,-1],
    crm.data['input_file']['rates']['recomb']['XS_cm2'][3,:,-1],
    label = 'bare->2P_1.5, ColRadPy'
    )

ax2.loglog(
    rr[1][1]['EGRID'],
    rr[1][1]['RR crosssection'][2,:]*1e-20,
    '*',
    label = 'bare->2P_1.5, FAC'
    )

ax2.set_title('H-like Ar, rad. recomb.')
ax2.grid('on')
leg = ax2.legend()
leg.set_draggable('on')
ax2.set_xlim(1e0, 1e5)
ax2.set_xlabel('E_incident [eV]')
ax2.set_ylabel('cross-section [cm2]')





ci = rfac.read_ci(file6+'a.ci')

E_grid = crm.data['input_file']['rates']['ioniz']['engy_eV'][0,:,-1]
# Electron kinetic momentum squared with fine structure correction
k02 = 2 * E_grid *eV2Hartree *(
    1+
    0.5 * fine_struct2 * E_grid *eV2Hartree
    ) # [atomic units]

# Formula for collision strength
xx1 = E_grid[:,None]/ci[1][0]['Delta E'][None,:]
xx2 = E_grid[:,None]/ci[1][1]['Delta E'][None,:]
yy1 = 1 - 1/xx1
yy2 = 1 - 1/xx2

omega1 = (
    ci[1][0]['parameters'][None,:,0] *np.log(xx1)
    + ci[1][0]['parameters'][None,:,1] *yy1**2
    + ci[1][0]['parameters'][None,:,2] *1/xx1 *yy1
    + ci[1][0]['parameters'][None,:,3] *1/xx1**2 *yy1
    )
omega2 = (
    ci[1][1]['parameters'][None,:,0] *np.log(xx2)
    + ci[1][1]['parameters'][None,:,1] *yy2**2
    + ci[1][1]['parameters'][None,:,2] *1/xx2 *yy2
    + ci[1][1]['parameters'][None,:,3] *1/xx2**2 *yy2
    )

# Cross-section
XS_tmp1 = (
    omega1
    / k02[:,None]
    / (1 +ci[1][0]['bound_2J'][None,:])
    ) * a02 # [cm2]
XS_tmp2 = (
    omega2
    / k02[:,None]
    / (1 +ci[1][1]['bound_2J'][None,:])
    ) * a02 # [cm2]



fig3, ax3 = plt.subplots()

ax3.loglog(
    E_grid,
    XS_tmp1[:,0],
    label = '(1s)2S_0.5->bare, formula'
    )

ax3.loglog(
    crm.data['input_file']['rates']['ioniz']['engy_eV'][0,:,-1],
    crm.data['input_file']['rates']['ioniz']['XS_cm2'][0,:,-1],
    label = '(1s)2S_0.5->bare, ColRadPy'
    )

ax3.loglog(
    ci[1][0]['EGRID']+ci[1][0]['TEGRID'],
    ci[1][0]['crosssection'][0,:]*1e-20,
    '*',
    label = '(1s)2S_0.5->bare, FAC'
    )

ax3.loglog(
    E_grid,
    XS_tmp2[:,0],
    label = '2P_0.5->bare, formula'
    )

ax3.loglog(
    crm.data['input_file']['rates']['ioniz']['engy_eV'][1,:,-1],
    crm.data['input_file']['rates']['ioniz']['XS_cm2'][1,:,-1],
    label = '2P_0.5->bare, ColRadPy'
    )

ax3.loglog(
    ci[1][1]['EGRID']+ci[1][1]['TEGRID'][0],
    ci[1][1]['crosssection'][0,:]*1e-20,
    '*',
    label = '2P_0.5->bare, FAC'
    )

ax3.loglog(
    E_grid,
    XS_tmp2[:,1],
    label = '(2s)2S_0.5->bare, formula'
    )

ax3.loglog(
    crm.data['input_file']['rates']['ioniz']['engy_eV'][2,:,-1],
    crm.data['input_file']['rates']['ioniz']['XS_cm2'][2,:,-1],
    label = '(2s)2S_0.5->bare, ColRadPy'
    )

ax3.loglog(
    ci[1][1]['EGRID']+ci[1][1]['TEGRID'][0],
    ci[1][1]['crosssection'][1,:]*1e-20,
    '*',
    label = '(2s)2S_0.5->bare, FAC'
    )

ax3.loglog(
    E_grid,
    XS_tmp2[:,2],
    label = '2P_1.5->bare, formula'
    )

ax3.loglog(
    crm.data['input_file']['rates']['ioniz']['engy_eV'][3,:,-1],
    crm.data['input_file']['rates']['ioniz']['XS_cm2'][3,:,-1],
    label = '2P_1.5->bare, ColRadPy'
    )

ax3.loglog(
    ci[1][1]['EGRID']+ci[1][1]['TEGRID'][0],
    ci[1][1]['crosssection'][2,:]*1e-20,
    '*',
    label = '2P_1.5->bare, FAC'
    )

ax3.set_title('H-like Ar, ioniz.')
ax3.grid('on')
leg = ax3.legend()
leg.set_draggable('on')
ax3.set_xlabel('E_incident [eV]')
ax3.set_ylabel('cross-section [cm2]')