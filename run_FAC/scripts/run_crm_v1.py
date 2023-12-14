'''

Test script for running collisional-radiatve modeling

cjperks

Dec. 8th, 2023

TO DO:
    1) Collisional ionization (rt6)
    2) Radiative recombination (rt4)
    3) Dielectronic capture (rt5)
    4) Higher energy levels beyond nmax=2 (cascades)
    5) Transport modeling in ionization balance equation
    6) Solve many charge states at once
    7) Couple to pfac.crm to get rates
    8_ Compare to Adam Foster's results


'''

# Modules
import aurora

# Plasma parameters
ne_cm3 = 1e14
Te_eV = 10e3
ni_cm3 = 1e10

ion = 'Ar'

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

# Path to FAC results
path_fac = '/home/cjperks/atomic_world/run_FAC/scripts/Ar/test'

# Charge states to include in modeling
csc = [17],

from pfac import rfac
import os

# Loads energy level data
en = rfac.read_en(
    os.path.join(path_fac, 'Ar01a.en')
    )

# Loads in rate data
rt1 = np.loadtxt(
    os.path.join(path_fac, 'Ar01a.rt1')
    ) # Spontaneous emissiion, dim(ntrans, 5), [1/s]

rt2 = np.loadtxt(
    os.path.join(path_fac, 'Ar01a.rt2')
    ) # 2-photon spontaneous emissiion, dim(ntrans, 5), [1/s]
if rt2.ndim == 1:
    rt2 = rt2[None,:]

rt3 = np.loadtxt(
    os.path.join(path_fac, 'Ar01a.rt3')
    ) # Collisional (de)excitation, dim(ntrans, 5), [1e-10 cm3/s]

rt4 = np.loadtxt(
    os.path.join(path_fac, 'Ar01a.rt4')
    ) # Radiative recombination, dim(ntrans, 5), [1e-10 cm3/s]

rt6 = np.loadtxt(
    os.path.join(path_fac, 'Ar01a.rt6')
    ) # Collisional ionization, dim(ntrans, 5), [1e-10 cm3/s]

##########################################
#
#           Organizes level data
#
##########################################

# Initializes dictionary to store level data
lvl = {}

# Maximum qunatum number to solve for energy levels
nmax = 2
lvls = []

# Loop over .en file blocks
for ind in np.arange(len(en[1])):

    # Loop over levels in block
    for ll, ilev in enumerate(en[1][ind]['ILEV']):
        # Dictionary to store data
        lvl[ilev] = {}

        # Level energy
        lvl[ilev]['E_eV'] = en[1][ind]['ENERGY'][ll]

        # Level qunatum numbers
        lvl[ilev]['VNL'] = en[1][ind]['VNL'][ll]

        # Level name
        lvl[ilev]['sname'] = en[1][ind]['sname'][ll]
        lvl[ilev]['name'] = en[1][ind]['name'][ll]

        # Number of electrons
        lvl[ilev]['nele'] = en[1][ind]['NELE']

        # If solving for this level population
        if ((int(lvl[ilev]['VNL']/100) <= nmax)
            and (lvl[ilev]['E_eV'] > 0)
            and (lvl[ilev]['nele'] == 1)):
            lvls.append(ilev)
            lvl[ilev]['ground'] = 'n'
        else:
            lvl[ilev]['ground'] = 'y'

        # --- Finds the spontaneous emission rates for this state --- #
        lvl[ilev]['spon_emis'] = {} # [1/s]

        # Loop over transitions in rt1
        for tran in np.arange(rt1.shape[0]):
            # If this is an upper state of the data
            if rt1[tran,1] == ilev:
                # stores data labeled by lower level
                lwr = rt1[tran,2]
                lvl[ilev]['spon_emis'][lwr] = rt1[tran,3]

        # --- Finds the 2-photon spontaneous emission rates for this state --- #

        # Loop over transitions in rt2
        for tran in np.arange(rt2.shape[0]):
            # If this is an upper state of the data
            if rt2[tran,1] == ilev:
                # stores data labeled by lower level
                lwr = rt2[tran,2]
                lvl[ilev]['spon_emis'][lwr] += rt2[tran,3]

        # --- Finds the collisional (de-)excitation for this state --- #
        lvl[ilev]['coll_excit'] = {} # [cm3/s]

        # Loop over transitions in rt3
        for tran in np.arange(rt3.shape[0]):
            # If this is a lower level of the data
            if rt3[tran,1] == ilev:
                # stores data labeled by upper level
                upr = rt3[tran,2]
                lvl[ilev]['coll_excit'][upr] = rt3[tran,3]*1e-10

        # --- Finds the collisional dexcitation for this state --- #
        lvl[ilev]['coll_de_excit'] = {} # [cm3/s]

        # Loop over transitions in rt3
        for tran in np.arange(rt3.shape[0]):
            # If this is a upper level of the data
            if rt3[tran,2] == ilev:
                # stores data labeled by upper level
                lwr = rt3[tran,1]
                lvl[ilev]['coll_de_excit'][lwr] = rt3[tran,4]*1e-10

        # --- Finds the radiative recombination for this state --- #
        lvl[ilev]['rad_recomb'] = {} # [cm3/s]

        # Loop over transitions in rt3
        for tran in np.arange(rt4.shape[0]):
            # If this is a bound level of the data
            if rt4[tran,2] == ilev:
                # stores data labeled by coninuum level
                cnt = rt4[tran,1]
                lvl[ilev]['rad_recomb'][cnt] = rt4[tran,3]*1e-10

            # If this is a continuum level of the data
            elif rt4[tran,1] == ilev:
                # stores data labeled by bound level
                bnd = rt4[tran,2]
                lvl[ilev]['rad_recomb'][bnd] = rt4[tran,3]*1e-10

        # --- Finds the collisional ionization for this state --- #
        lvl[ilev]['coll_ion'] = {} # [cm3/s]

        # Loop over transitions in rt3
        for tran in np.arange(rt6.shape[0]):
            # If this is a bound level of the data
            if rt6[tran,1] == ilev:
                # stores data labeled by coninuum level
                cnt = rt6[tran,2]
                lvl[ilev]['coll_ion'][cnt] = rt6[tran,3]*1e-10
        


##########################################
#
#      Rate eqn matrix
#
##########################################

# Initializes matrix
mat = np.zeros((len(lvls), len(lvls))) # dim(nlvls, nlvls)
rhs = np.zeros(len(lvls))

def _calc_cascade(
    quant = 'coll_excit',
    lvl = None,
    # Levels in consideration
    lwr = None,
    upr = None,
    # Plasma
    ne_cm3 = None,
    nmax = None,
    ):

    data = lvl[lwr][quant]

    # Initializes cascade term
    casc = 0

    # Loop over lvls
    for higher in data.keys():
        # If higher energy state
        if ((lvl[upr]['E_eV'] < lvl[higher]['E_eV'])
            and (int(lvl[higher]['VNL']/100) <= nmax)):
            # Gets data from lower to higher
            data_hu = data[higher]

            # Calculates branching ratio from higher to upper
            b_hu_top = (
                lvl[higher]['spon_emis'][upr]
                + ne_cm3 * lvl[higher]['coll_de_excit'][upr]
                )
            b_hu_bot  = 0
            for key in lvl[higher]['spon_emis'].keys():
                b_hu_bot += lvl[higher]['spon_emis'][key]
            for key in lvl[higher]['coll_excit'].keys():
                b_hu_bot += ne_cm3 *lvl[higher]['coll_excit'][key]
            for key in lvl[higher]['coll_de_excit'].keys():
                b_hu_bot += ne_cm3 *lvl[higher]['coll_de_excit'][key]

            casc += data_hu *b_hu_top/b_hu_bot

    return casc *ne_cm3

'''

Population equation into state rr of charge Z

Assume:
    1) Optical thin 
        -> Ignore radiation field interactions
    2) Lower energy density 
        -> Dominant population is ground state for all charge states
    3) Quasi-equilibrium
        -> Transport time scales and atomic time scales separable
        -> Can solve ionization balance + transport equation first with Aurora
        -> dot{N_rr^Z} = 0

Solve:
dot{N_rr^Z} =  (
    N_rr^Z * sink                           # Sink term
    + \sum_{E_ll < E_rr} N_ll^Z * src_lwr   # Source term, 
                                                # from lower energy states
    + \sum_{E_ll > E_rr} N_ll^Z * src_upr   # Source term, 
                                                # from higher energy states
    + src_gnd                               # Source term,
                                                # exciting charge Z ground
    + src_ion                               # Source term,
                                                # ionizing charge Z-1
    + src_rec                               # Source term,
                                                # recombining charge Z+1
    
sink = -1*(
    \sum_{E_ll < E_rr} A_{rr->ll}       # Decay to lower states
    + N_e * S_{rr -> cnt}               # Collisional ionization to continuum
    + N_e * \sum{ll=/=rr} S_{rr->ll}    # Collisional excitation
    )                                   # [1/s]

src_lwr = N_e * S_{ll->rr}^             # Collisional excitation (excl. ground)
                                        # Incl. Cascasde from E_kk > E_rr
                                        # [1/s]

src_upr = (
    N_e * S_{ll->rr}                    # Collisional excitation
                                        # [1/s]
    + A_{ll->rr}                        # Spontaneous emission
    )

src_gnd = N_e * N_gg^Z * S_{ll->rr}^    # Collisional excitation from ground
                                        # Incl. Cascasde from E_kk > E_rr
                                        # [1/cm3/s]

src_ion = N_e * N_gg^{Z-1} * S_{gg->rr} # Collisional ionization from Z-1 ground
                                        # [1/cm3/s]

src_rec = (
    N_e * N_gg^{Z+1} * ar_{gg->rr}^     # Radiative recombination from Z+1 ground
                                        # Incl. Cascade from E_kk > E_rr
    + N_e * N_gg^{Z+1} * ad_{gg->rr}^   # Dielectronic capture from Z+1 ground
                                        # Incl. Cascade from E_kk > E_rr
    )                                   # [1/cm3/s]

'''

# Calculates the sink term in population equation
def _calc_sink(
    lvl = None,
    state = None,
    ne_cm3 = None,
    ):

    # Init
    sink = 0.0

    # Radiative decay out of state rr
    for key in lvl[state]['spon_emis'].keys():
        sink -= lvl[state]['spon_emis'][key]

    # Collisional excitation/de-excitation out of state rr
    for key in lvl[state]['coll_excit'].keys():
        sink -= ne_cm3 *lvl[state]['coll_excit'][key]
    for key in lvl[state]['coll_de_excit'].keys():
        sink -= ne_cm3 *lvl[state]['coll_de_excit'][key]

    # Collisional ionization out of state rr
    for key in lvl[state]['coll_ion'].keys():
        sink -= ne_cm3 *lvl[lvls[rr]]['coll_ion'][key]

    # [1/s]
    return sink

# Source from lower energy states
def _calc_src_lwr(
    lvl = None,
    lwr = None,
    upr = None,
    ne_cm3 = None,
    nmax = None,
    ):

    # Init
    src = 0.0

    # Collisional excitation from cc into rr including cascades
    src += ne_cm3 *lvl[lwr]['coll_excit'][upr]

    src += _calc_cascade(
        quant = 'coll_excit',
        lvl = lvl,
        lwr = lwr,
        upr = upr,
        ne_cm3=ne_cm3,
        nmax=nmax,
        )

    # [1/s]
    return src

# Source from higher energy states
def _calc_src_upr(
    lvl = None,
    upr = None,
    lwr = None,
    ne_cm3 = None,
    ):

    # Init
    src = 0.0

    # Collisional de-excitation from cc to rr
    src += ne_cm3 *lvl[upr]['coll_de_excit'][lwr]

    # Radiative decay from cc to rr
    src += lvl[upr]['spon_emis'][lwr]

    # [1/s]
    return src

# Source from charge Z ground state
def _calc_src_gnd(
    lvl = None,
    gnd = None,
    upr = None,
    ne_cm3 = None,
    nz_cm3 = None,
    ZZ = None,
    nmax = None,
    ):

    # Init
    src = 0.0

    # Ground state of ion Z collisional excite to rr including cascades
    src += ne_cm3 *lvl[gnd]['coll_excit'][upr] *nz_cm3[ZZ]

    src += nz_cm3[17] * _calc_cascade(
        quant = 'coll_excit',
        lvl = lvl,
        lwr = gnd,
        upr = upr,
        ne_cm3 = ne_cm3,
        nmax = nmax,
        )

    # [1/cm3/s]
    return src

# Source from charge Z+1 ground state recombining to state rr
def _calc_src_rec(
    lvl = None,
    state = None,
    ne_cm3 = None,
    nz_cm3 = None,
    ZZ = None,
    nmax = None,
    ):
    '''
    TO DO:
        1) Add dielectronic recombination
    '''

    # Init 
    src = 0.0

    # Ground state of ion Z+1 radiative recombines to rr including cascades
    for key in lvl[state]['rad_recomb'].keys():
        src += ne_cm3 *lvl[state]['rad_recomb'][key] * nz_cm3[ZZ +1]

        src += _calc_cascade(
                quant = 'rad_recomb',
                lvl = lvl,
                lwr = lvl[key],
                upr = state,
                ne_cm3=ne_cm3,
                nmax=nmax,
                ) *nz_cm3[ZZ +1]


    # [1/cm3/s]
    return src

# Source from charge Z+1 ground state ionizing to state rr
def _calc_src_ion(
    lvl = None,
    ):
    '''
    TO DO:
        1) Add collisional ionization
    '''

    # Init 
    src = 0.0

    # [1/cm3/s]
    return src


# Loop over matrix rows -> the specific rate equation for state rr
for rr in np.arange(mat.shape[0]):
    # Loop over matrix columns -> Sources from state cc
    for cc in np.arange(mat.shape[1]):
        # If considering the same state
        if rr == cc:
            mat[rr,cc] += _calc_sink(
                lvl=lvl,
                state=lvls[rr],
                ne_cm3=ne_cm3,
                ) # [1/s]

        # If considering a lower energy state 
        elif lvl[lvls[cc]]['E_eV'] < lvl[lvls[rr]]['E_eV']:
            mat[rr,cc] += _calc_src_lwr(
                lvl=lvl,
                lwr=lvls[cc],
                upr=lvls[rr],
                ne_cm3=ne_cm3,
                nmax=nmax,
                )

        # If considering a higher energy state
        elif lvl[lvls[cc]]['E_eV'] > lvl[lvls[rr]]['E_eV']:
            mat[rr,cc] += _calc_src_lwr(
                lvl=lvl,
                upr=lvls[cc],
                lwr=lvls[rr],
                ne_cm3=ne_cm3,
                )

    # Ground state of ion Z collisional excite to rr including cascades
    rhs[rr] += _calc_src_gnd(
        lvl = lvl,
        gnd = 0,
        upr = lvls[rr],
        ne_cm3 = ne_cm3,
        nz_cm3 = nz_cm3,
        ZZ = ZZ,
        nmax = namx,
        ) # [1/cm3/s]

    # Ground state of ion Z+1 radiative recombines to rr including cascades
    rhs[rr] += _calc_src_rec(
        lvl = lvl,
        state = lvls[rr],
        ne_cm3 = ne_cm3,
        nz_cm3 = nz_cm3,
        ZZ = ZZ,
        nmax = nmax,
        ) # [1/cm3/s]

    # Ground state of ion Z-1 collisionally ionizes to rr including cascades
    rhs[rr] += _calc_src_ion(
        lvl = lvl
        ) # [1/cm3/s]

# Solves for the population densities
pops = np.linalg.solve(mat, -1*rhs) # [1/cm3]
