'''

Debugging differences in excitation rate coefficients converted
to upsilon form

cjperks
Feb 28, 2024

'''

# Modules
import sys, os
from colradpy import colradpy

# Ion
sp = 'Ar'
nele = 1
Znuc = 18

# Common FAC atomic data files name
file_me = os.path.join(
    '/home/cjperks/2008_SPARCxray',
    'FAC/ATOM',
    sp,
    sp+'%02d'%(nele)
    )

file_04 = os.path.join(
    '/home/cjperks',
    '2008_SPARCxray/ColRadPy',
    'scripts/learning',
    'ssh42_cs_ic#ar17.dat'
    )

# Simulation grids
temp_arr = np.logspace(np.log10(8.62), np.log10(8.5e4), 41) # [eV]
dens_arr = np.logspace(np.log10(1e10), np.log10(1e20), 21) # [cm^-3]
meta_arr = np.array([0])

# Load atomic data
me = colradpy(
    file_me,
    meta_arr,
    temp_arr,
    dens_arr,
    use_recombination=False,
    use_recombination_three_body = False,
    use_ionization=False,
    suppliment_with_ecip=False,
    # New material
    atomic_data_type = 'FAC',
    ele = sp,
    nele = nele,
    Zele = Znuc,
    EEDF = 'Maxwellian',
    atomic_physics = 'incl_all',
    )

# Load atomic data
adf = colradpy(
    file_04,
    meta_arr,
    temp_arr,
    dens_arr,
    use_recombination=False,
    use_recombination_three_body = False,
    use_ionization=False,
    suppliment_with_ecip=False,
    )


ind_adf = 5
ind_me = 5

tr_me = me.data['rates']['excit']['col_transitions'][ind_me,:]
tr_adf = adf.data['rates']['excit']['col_transitions'][ind_adf,:]

print('me')
print(tr_me)
print('upr')
print(me.data['atomic']['S'][tr_me[0]-1])
print(me.data['atomic']['L'][tr_me[0]-1])
print(me.data['atomic']['w'][tr_me[0]-1])
print('lwr')
print(me.data['atomic']['S'][tr_me[1]-1])
print(me.data['atomic']['L'][tr_me[1]-1])
print(me.data['atomic']['w'][tr_me[1]-1])

print('adf04')
print(tr_adf)
print('upr')
print(adf.data['atomic']['S'][tr_adf[0]-1])
print(adf.data['atomic']['L'][tr_adf[0]-1])
print(adf.data['atomic']['w'][tr_adf[0]-1])
print('lwr')
print(adf.data['atomic']['S'][tr_adf[1]-1])
print(adf.data['atomic']['L'][tr_adf[1]-1])
print(adf.data['atomic']['w'][tr_adf[1]-1])

fig, ax = plt.subplots()

eV2K = 11604.5

ax.loglog(
    me.data['input_file']['temp_grid']/eV2K,
    me.data['rates']['excit']['col_excit'][ind_me,:],
    label = 'me'
    )

ax.loglog(
    adf.data['input_file']['temp_grid']/eV2K,
    adf.data['rates']['excit']['col_excit'][ind_adf,:],
    label = 'adf04'
    )

ax.legend()

ind_me = 7

plt.loglog(
    me.data['rates']['excit']['engy_eV'][ind_me,:,-1],
    me.data['rates']['excit']['XS_cm2'][ind_me,:,-1]
    )