'''

Script to validate modeling of He-like Ar

1) the 1s2-1s2l transitions and
2) the 1snl-2l'nl satellite transitions

March 6th, 2024


'''



# Module
import os, sys
import aurora

# Ion
sp = 'Ar'
nele = 2
Znuc = 18
ion_A = 40
fi = 1e-3

# Plasma
Te_eV = 2e3
Ti_eV = 2e3
ne_cm3 = 1e14

# Files
file_me = os.path.join(
    '/home/cjperks',
    '2008_SPARCxray/ColRadPy/FAC',
    'test_FAC_'+sp+'%02d.dat'%(nele)
    )

file_AF = os.path.join(
    '/home/cjperks',
    'tofu_sparc/atomic_data',
    'ADAS_PEC_files/fsciortino_share/atomdb',
    'pec#'+sp.lower()+'%02d.dat'%(Znuc-nele)
    )

# Broadening mechanisms
dbroad_AF_G = {
    'Voigt':{
        'Ti_eV': Ti_eV,
        'ion_A': ion_A,
        'key_options': 'wavelength',
        }
    }
dbroad_AF_V = {
    'Voigt':{
        'Ti_eV': Ti_eV,
        'ion_A': ion_A,
        'key_options': 'wavelength',
        '3.994153': 4.431127e6,       # [1/s], z
        '3.969363': 1.637323e12,      # [1/s], y
        '3.965866': 3.086273e8,       # [1/s], x
        '3.949075': 1.097426e14,      # [1/s], w
        }
    }

dbroad_me= {
    'Voigt':{
        'Ti_eV': Ti_eV,
        'ion_A': ion_A,
        'key_options': 'wavelength',
        '3.99371': 4.431127e6,       # [1/s], z
        '3.96930': 1.637323e12,      # [1/s], y
        '3.96582': 3.086273e8,       # [1/s], x
        '3.94901': 1.097426e14,      # [1/s], w
        },
    '2-photon':{
        'wavelength': [3.96711],
        'Znuc': 18,
        'nele': 2,
        },
    }



# Charge states per ionization balance
atom_data = aurora.atomic.get_atom_data(sp,['scd','acd'])
_, fz = aurora.atomic.get_frac_abundances(
    atom_data, 
    ne_cm3, Te_eV, 
    plot=False,
    )

# Aurora spectrum modeling
out_AF_G = aurora.get_local_spectrum(
    file_AF,
    ne_cm3,
    Te_eV,
    ion_exc_rec_dens = [
        fi*ne_cm3*fz[0,Znuc-nele-1],
        fi*ne_cm3*fz[0,Znuc-nele],
        fi*ne_cm3*fz[0,Znuc-nele+1]
        ],
    Ti_eV=Ti_eV,
    n0_cm3 = 0.0,
    dlam_A=0.0,
    plot_spec_tot=False,
    # New Material
    dbroad=dbroad_AF_G
    )

out_AF_V = aurora.get_local_spectrum(
    file_AF,
    ne_cm3,
    Te_eV,
    ion_exc_rec_dens = [
        fi*ne_cm3*fz[0,Znuc-nele-1],
        fi*ne_cm3*fz[0,Znuc-nele],
        fi*ne_cm3*fz[0,Znuc-nele+1]
        ],
    Ti_eV=Ti_eV,
    n0_cm3 = 0.0,
    dlam_A=0.0,
    plot_spec_tot=False,
    # New Material
    dbroad=dbroad_AF_V
    )

out_me = aurora.get_local_spectrum(
    file_me,
    ne_cm3,
    Te_eV,
    ion_exc_rec_dens = [
        fi*ne_cm3*fz[0,Znuc-nele-1],
        fi*ne_cm3*fz[0,Znuc-nele],
        fi*ne_cm3*fz[0,Znuc-nele+1]
        ],
    Ti_eV=Ti_eV,
    n0_cm3 = 0.0,
    dlam_A=0.0,
    plot_spec_tot=False,
    # New Material
    dbroad=dbroad_me
    )




fig,ax = plt.subplots(1,2)

ax[0].plot(
    out_AF_G[0],
    out_AF_G[1],
    label = 'AF, ion, Gauss'
    )
ax[0].plot(
    out_AF_G[0],
    out_AF_G[2],
    label = 'AF, exc, Gauss'
    )
ax[0].plot(
    out_AF_G[0],
    out_AF_G[3],
    label = 'AF, RR, Gauss'
    )

ax[0].plot(
    out_AF_V[0],
    out_AF_V[1],
    label = 'AF, ion, Voigt'
    )
ax[0].plot(
    out_AF_V[0],
    out_AF_V[2],
    label = 'AF, exc, Voigt'
    )
ax[0].plot(
    out_AF_V[0],
    out_AF_V[3],
    label = 'AF, RR, Voigt'
    )

ax[0].plot(
    out_me[0],
    out_me[2],
    label = 'me, exc, Voigt'
    )

if out_me[3] is not None:
    ax[0].plot(
        out_me[0],
        out_me[3],
        label = 'me, rec, Voigt'
        )

ax[0].set_xlim(3.94, 4.00)
leg = ax[0].legend()
leg.set_draggable('on')
ax[0].grid('on')
ax[0].set_ylabel('emissivity [ph/s/cm3/AA]')
ax[0].set_xlabel('wavelength [AA]')
ax[0].set_title(sp+'%02d+, Te = %01.1f keV, Ti = %01.1f keV'%(Znuc-nele, Te_eV/1e3, Ti_eV/1e3))



ax[1].plot(
    out_AF_G[0],
    out_AF_G[1],
    label = 'AF, ion, Gauss'
    )
ax[1].plot(
    out_AF_G[0],
    out_AF_G[2],
    label = 'AF, exc, Gauss'
    )
ax[1].plot(
    out_AF_G[0],
    out_AF_G[3],
    label = 'AF, RR, Gauss'
    )
ax[1].plot(
    out_AF_G[0],
    out_AF_G[4],
    label = 'AF, DR, Gauss'
    )


ax[1].plot(
    out_AF_V[0],
    out_AF_V[1],
    label = 'AF, ion, Voigt'
    )
ax[1].plot(
    out_AF_V[0],
    out_AF_V[2],
    label = 'AF, exc, Voigt'
    )
ax[1].plot(
    out_AF_V[0],
    out_AF_V[3],
    label = 'AF, RR, Voigt'
    )

ax[1].plot(
    out_AF_V[0],
    out_AF_V[4],
    label = 'AF, DR, Voigt'
    )

ax[1].plot(
    out_me[0],
    out_me[2],
    label = 'me, exc, Voigt'
    )

if out_me[3] is not None:
    ax[1].plot(
        out_me[0],
        out_me[3],
        label = 'me, rec, Voigt'
        )

ax[1].set_xlim(3.725, 3.745)
leg = ax[1].legend()
leg.set_draggable('on')
ax[1].grid('on')
ax[1].set_ylabel('emissivity [ph/s/cm3/AA]')
ax[1].set_xlabel('wavelength [AA]')
ax[1].set_title(sp+'%02d+, Te = %01.1f keV, Ti = %01.1f keV'%(Znuc-nele, Te_eV/1e3, Ti_eV/1e3))
ax[1].set_ylim(0,5e14)


##################################################################


# Modules
import sys, os
from colradpy import colradpy

# Ion
sp = 'Ar'
nele = 2
Znuc = 18

# Common FAC atomic data files name
file = os.path.join(
    '/home/cjperks/2008_SPARCxray',
    #'FAC/ATOM',
    #sp,
    'FAC/ATOM_editted_v1',
    sp+'_v11',
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
'''
ind_c = 5-5
ind_f = 5-5
plt.loglog(
    crm.data['input_file']['rates']['excit']['engy_eV'][ind_c,:,-1],
    crm.data['input_file']['rates']['excit']['XS_cm2'][ind_c,:,-1]
    )
plt.loglog(
    ce[1][0]['EGRID']+ce[1][0]['TE0'],
    ce[1][0]['crosssection'][ind_f,:]*1e-20,
    '*'
    )
'''
################################################################3
temp_arr = np.logspace(np.log10(8.62), np.log10(8.5e4), 41) # [eV]


# Reads ADF15
trs_AF = aurora.read_adf15(file_AF)
trs_me = aurora.read_adf15(file_me)

fig1,ax1 = plt.subplots(2,2)
fig1.tight_layout(pad=1.5)


# Line of interest
lamb_AF = 3.949075
inds = np.where(trs_AF['lambda [A]'] == lamb_AF)

ind_AF_i = inds[0][0]
ind_AF_e = inds[0][1]
ind_AF_r = inds[0][2]

# PECs
val_AF_i = 10**trs_AF['log10 PEC fun'][ind_AF_i].ev(np.log10(1e14), np.log10(temp_arr))
val_AF_e = 10**trs_AF['log10 PEC fun'][ind_AF_e].ev(np.log10(1e14), np.log10(temp_arr))
val_AF_r = 10**trs_AF['log10 PEC fun'][ind_AF_r].ev(np.log10(1e14), np.log10(temp_arr))



#lamb_me = 3.73612
lamb_me = 3.94901
inds = np.where(trs_me['lambda [A]'] == lamb_me)

ind_me_e = inds[0][0]
ind_me_r = inds[0][1]
#ind_me_r = None


val_me_ew = 10**trs_me['log10 PEC fun'][ind_me_e].ev(np.log10(1e14), np.log10(temp_arr))
if ind_me_r is not None:
    val_me_rw = 10**trs_me['log10 PEC fun'][ind_me_r].ev(np.log10(1e14), np.log10(temp_arr))




ax1[0,0].plot(
    temp_arr/1e3,
    val_AF_i,
    label = 'AF, ioniz'
    )
ax1[0,0].plot(
    temp_arr/1e3,
    val_AF_e,
    label = 'AF, excit'
    )
ax1[0,0].plot(
    temp_arr/1e3,
    val_AF_r,
    label = 'AF, rad recomb'
    )


ax1[0,0].plot(
    temp_arr/1e3,
    val_me_ew,
    label = 'me, excit'
    )
if ind_me_r is not None:
    ax1[0,0].plot(
        temp_arr/1e3,
        val_me_rw,
        label = 'me, recomb'
        )
'''
ax1[0,0].plot(
    temp_arr/1e3,
    crm.data['input_file']['rates']['excit']['ratec_cm3/s'][2,:],
    label = 'me, excit, rate'
    )

ax1[0,0].plot(
    temp_arr/1e3,
    crm.data['input_file']['rates']['recomb']['recomb_excit'][3,:],
    label = 'me, recomb, rate'
    )
'''
ax1[0,0].grid('on')
ax1[0,0].set_xlabel(r'$T_e$ [$keV$]')
ax1[0,0].set_ylabel(r'PEC [$cm^3/s$]')
ax1[0,0].set_title('w')
leg = ax1[0,0].legend()
leg.set_draggable('on')
#ax1[0,0].set_ylim(0,5e-12)




# Line of interest
lamb_AF = 3.965866
inds = np.where(trs_AF['lambda [A]'] == lamb_AF)

ind_AF_i = inds[0][0]
ind_AF_e = inds[0][1]
ind_AF_r = inds[0][2]


# PECs
val_AF_i = 10**trs_AF['log10 PEC fun'][ind_AF_i].ev(np.log10(1e14), np.log10(temp_arr))
val_AF_e = 10**trs_AF['log10 PEC fun'][ind_AF_e].ev(np.log10(1e14), np.log10(temp_arr))
val_AF_r = 10**trs_AF['log10 PEC fun'][ind_AF_r].ev(np.log10(1e14), np.log10(temp_arr))



#lamb_me = 3.73612
lamb_me = 3.96582
inds = np.where(trs_me['lambda [A]'] == lamb_me)

ind_me_e = inds[0][0]
ind_me_r = inds[0][1]
#ind_me_r = None


val_me_ex = 10**trs_me['log10 PEC fun'][ind_me_e].ev(np.log10(1e14), np.log10(temp_arr))
if ind_me_r is not None:
    val_me_rx = 10**trs_me['log10 PEC fun'][ind_me_r].ev(np.log10(1e14), np.log10(temp_arr))



ax1[0,1].plot(
    temp_arr/1e3,
    val_AF_i,
    label = 'AF, ioniz'
    )
ax1[0,1].plot(
    temp_arr/1e3,
    val_AF_e,
    label = 'AF, excit'
    )
ax1[0,1].plot(
    temp_arr/1e3,
    val_AF_r,
    label = 'AF, rad recomb'
    )


ax1[0,1].plot(
    temp_arr/1e3,
    val_me_ex,
    label = 'me, excit'
    )
if ind_me_r is not None:
    ax1[0,1].plot(
        temp_arr/1e3,
        val_me_rx,
        label = 'me, recomb'
        )
'''
ax1[0,0].plot(
    temp_arr/1e3,
    crm.data['input_file']['rates']['excit']['ratec_cm3/s'][2,:],
    label = 'me, excit, rate'
    )

ax1[0,0].plot(
    temp_arr/1e3,
    crm.data['input_file']['rates']['recomb']['recomb_excit'][3,:],
    label = 'me, recomb, rate'
    )
'''
ax1[0,1].grid('on')
ax1[0,1].set_xlabel(r'$T_e$ [$keV$]')
ax1[0,1].set_ylabel(r'PEC [$cm^3/s$]')
ax1[0,1].set_title('x')
leg = ax1[0,1].legend()
leg.set_draggable('on')
#ax1[0,1].set_ylim(0,5e-12)



# Line of interest
lamb_AF = 3.969363
inds = np.where(trs_AF['lambda [A]'] == lamb_AF)

ind_AF_i = inds[0][0]
ind_AF_e = inds[0][1]
ind_AF_r = inds[0][2]

# PECs
val_AF_i = 10**trs_AF['log10 PEC fun'][ind_AF_i].ev(np.log10(1e14), np.log10(temp_arr))
val_AF_e = 10**trs_AF['log10 PEC fun'][ind_AF_e].ev(np.log10(1e14), np.log10(temp_arr))
val_AF_r = 10**trs_AF['log10 PEC fun'][ind_AF_r].ev(np.log10(1e14), np.log10(temp_arr))


#lamb_me = 3.73612
lamb_me = 3.96930
inds = np.where(trs_me['lambda [A]'] == lamb_me)

ind_me_e = inds[0][0]
ind_me_r = inds[0][1]
#ind_me_r = None


val_me_ey = 10**trs_me['log10 PEC fun'][ind_me_e].ev(np.log10(1e14), np.log10(temp_arr))
if ind_me_r is not None:
    val_me_ry = 10**trs_me['log10 PEC fun'][ind_me_r].ev(np.log10(1e14), np.log10(temp_arr))




ax1[1,0].plot(
    temp_arr/1e3,
    val_AF_i,
    label = 'AF, ioniz'
    )
ax1[1,0].plot(
    temp_arr/1e3,
    val_AF_e,
    label = 'AF, excit'
    )
ax1[1,0].plot(
    temp_arr/1e3,
    val_AF_r,
    label = 'AF, rad recomb'
    )


ax1[1,0].plot(
    temp_arr/1e3,
    val_me_ey,
    label = 'me, excit'
    )
if ind_me_r is not None:
    ax1[1,0].plot(
        temp_arr/1e3,
        val_me_ry,
        label = 'me, recomb'
        )
'''
ax1[1,0].plot(
    temp_arr/1e3,
    crm.data['input_file']['rates']['excit']['ratec_cm3/s'][2,:],
    label = 'me, excit, rate'
    )

ax1[1,0].plot(
    temp_arr/1e3,
    crm.data['input_file']['rates']['recomb']['recomb_excit'][3,:],
    label = 'me, recomb, rate'
    )
'''
ax1[1,0].grid('on')
ax1[1,0].set_xlabel(r'$T_e$ [$keV$]')
ax1[1,0].set_ylabel(r'PEC [$cm^3/s$]')
ax1[1,0].set_title('y')
leg = ax1[1,0].legend()
leg.set_draggable('on')
#ax1[1,0].set_ylim(0,5e-12)





# Line of interest
lamb_AF = 3.994153
inds = np.where(trs_AF['lambda [A]'] == lamb_AF)

ind_AF_i = inds[0][0]
ind_AF_e = inds[0][1]
ind_AF_r = inds[0][2]

# PECs
val_AF_i = 10**trs_AF['log10 PEC fun'][ind_AF_i].ev(np.log10(1e14), np.log10(temp_arr))
val_AF_e = 10**trs_AF['log10 PEC fun'][ind_AF_e].ev(np.log10(1e14), np.log10(temp_arr))
val_AF_r = 10**trs_AF['log10 PEC fun'][ind_AF_r].ev(np.log10(1e14), np.log10(temp_arr))


#lamb_me = 3.73612
lamb_me = 3.99371
inds = np.where(trs_me['lambda [A]'] == lamb_me)

ind_me_e = inds[0][0]
ind_me_r = inds[0][1]
#ind_me_r = None


val_me_ez = 10**trs_me['log10 PEC fun'][ind_me_e].ev(np.log10(1e14), np.log10(temp_arr))
if ind_me_r is not None:
    val_me_rz = 10**trs_me['log10 PEC fun'][ind_me_r].ev(np.log10(1e14), np.log10(temp_arr))




ax1[1,1].plot(
    temp_arr/1e3,
    val_AF_i,
    label = 'AF, ioniz'
    )
ax1[1,1].plot(
    temp_arr/1e3,
    val_AF_e,
    label = 'AF, excit'
    )
ax1[1,1].plot(
    temp_arr/1e3,
    val_AF_r,
    label = 'AF, rad recomb'
    )


ax1[1,1].plot(
    temp_arr/1e3,
    val_me_ez,
    label = 'me, excit'
    )
if ind_me_r is not None:
    ax1[1,1].plot(
        temp_arr/1e3,
        val_me_rz,
        label = 'me, recomb'
        )
'''
ax1[1,1].plot(
    temp_arr/1e3,
    crm.data['input_file']['rates']['excit']['ratec_cm3/s'][2,:],
    label = 'me, excit, rate'
    )

ax1[1,1].plot(
    temp_arr/1e3,
    crm.data['input_file']['rates']['recomb']['recomb_excit'][3,:],
    label = 'me, recomb, rate'
    )
'''
ax1[1,1].grid('on')
ax1[1,1].set_xlabel(r'$T_e$ [$keV$]')
ax1[1,1].set_ylabel(r'PEC [$cm^3/s$]')
ax1[1,1].set_title('z')
leg = ax1[1,1].legend()
leg.set_draggable('on')
#ax1[1,0].set_ylim(0,5e-12)

'''
ax1[0,0].set_xscale('log')
ax1[0,1].set_xscale('log')
ax1[1,0].set_xscale('log')
ax1[1,1].set_xscale('log')

ax1[0,0].set_yscale('log')
ax1[0,1].set_yscale('log')
ax1[1,0].set_yscale('log')
ax1[1,1].set_yscale('log')

ax1[0,0].set_ylim(1e-20,1e-10)
ax1[0,1].set_ylim(1e-20,1e-10)
ax1[1,0].set_ylim(1e-20,1e-10)
ax1[1,1].set_ylim(1e-20,1e-10)
'''