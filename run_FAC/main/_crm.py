'''

_crm.py is a module meant to preform collisional radiative modeling
in an optically thin plasma using Aurora and FAC modeling. 

Assumptions:
    1) Transport modeling already done with Aurora
    2) Any population with n>n_max is zero
    3) Radiation field negligible

TO DO:
    1) Write PEC files -> loop over ne, Te
    2) Compare against Adam's files
    3) Add more states in the data files
    4) More efficient way to get Maxwell rates

cjperks

Dec. 8th, 2023


'''

# Modules
from pfac import rfac
import os
import numpy as np
from scipy.interpolate import interp1d

__all__ = [
    'prep',
    'pops',
    'PECs',
    ]

#########################################
#
#           Main
#
#########################################

# Performs population balance
def pops(
    data = None,
    lvls = None,    # [list], levels to solve for
    nsim = None,    # Maximum qunatum number to incl. cascades
    # Plasma
    ne_cm3 = None,
    nz_cm3 = None,  # dim(cs,)
    ZZ = None,
    ):

    # Initializes matrix
    mat = np.zeros((len(lvls), len(lvls))) # dim(nlvls, nlvls)
    rhs = np.zeros(len(lvls))

    # Loop over matrix rows -> the specific rate equation for state rr
    for rr in np.arange(mat.shape[0]):
        # Loop over matrix columns -> Sources from state cc
        for cc in np.arange(mat.shape[1]):
            # If considering the same state
            if rr == cc:
                mat[rr,cc] += _calc_sink(
                    data=data,
                    state=lvls[rr],
                    ne_cm3=ne_cm3,
                    ) # [1/s]

            # If considering a lower energy state, including cascades
            elif data[lvls[cc]]['E_eV'] < data[lvls[rr]]['E_eV']:
                mat[rr,cc] += _calc_src_lwr(
                    data=data,
                    lwr=lvls[cc],
                    upr=lvls[rr],
                    ne_cm3=ne_cm3,
                    nsim=nsim,
                    ) # [1/s]

            # If considering a higher energy state
            elif data[lvls[cc]]['E_eV'] > data[lvls[rr]]['E_eV']:
                mat[rr,cc] += _calc_src_upr(
                    data=data,
                    upr=lvls[cc],
                    lwr=lvls[rr],
                    ne_cm3=ne_cm3,
                    )

        # Ground state of ion Z+1 radiative recombines to rr including cascades
        rhs[rr] += _calc_src_rec(
            data = data,
            state = lvls[rr],
            ne_cm3 = ne_cm3,
            nz_cm3 = nz_cm3,
            ZZ = ZZ,
            nsim = nsim,
            ) # [1/cm3/s]

        # Ground state of ion Z-1 collisionally ionizes to rr including cascades
        rhs[rr] += _calc_src_ion(
            data = data,
            ) # [1/cm3/s]

    # Solves for the population densities
    pops = np.linalg.solve(mat, -1*rhs) # [1/cm3]

    return pops/pops[0], mat, rhs


# Determines PEC terms given soln of population balance
def PECs(
    # States for transition of interest
    state = None,
    final = 0,
    # Data
    data = None,
    pops = None,
    nsim = None,    # Maximum qunatum number to incl. cascades
    # Plasma
    ne_cm3 = None,
    nz_cm3 = None,
    ZZ = None,
    ):

    # Excitation, [cm3/s]
    exc = pops[state]*data[state]['spon_emis'][final]/ne_cm3
 
    # Recombination from Z+1 ground
    rec = _calc_src_rec(
        data = data,
        state = state,
        ne_cm3 = ne_cm3,
        nz_cm3 = nz_cm3,
        ZZ = ZZ,
        nsim = nsim,
        ) # [1/cm3/s]

    # Ionization from Z-1 ground
    ion = _calc_src_ion(
        data = data,
        ) # [1/cm3/s]

    # Sink term
    sink =  -1*_calc_sink(
        data=data,
        state=state,
        ne_cm3=ne_cm3,
        ) # [1/s]


    # Intensity, [cm3/s]
    rec /= nz_cm3[ZZ+1]*ne_cm3
    ion /= nz_cm3[ZZ-1]*ne_cm3
    rec *= data[state]['spon_emis'][final]/sink
    ion *= data[state]['spon_emis'][final]/sink

    # Output, [cm3/s]
    return exc, rec, ion

#########################################
#
#           Data Management
#
#########################################

# Format FAC data
def prep(
    path_fac = None,
    ion = None,
    nele = None,
    nsim = None, # Maximum quantum number to solve for energy levels
    mode = None,
    Te_eV = None,
    ):

    # Data files
    if mode == 'rt':
        (
            en, rt1, rt2,
            rt3, rt4, rt6
            ) = _read_FAC_rt(
                path_fac = path_fac,
                ion = ion,
                nele = nele,
                )
    elif mode == 'mr':
        (
            en, tr, ce,
            rr, ci
            ) = _read_FAC_mr(
                path_fac = path_fac,
                ion = ion,
                nele = nele,
                )

    # Initializes dictionary to store level data
    data = {}
    lvls = []

    # Loop over .en file blocks
    for ind in np.arange(len(en[1])):

        # Loop over levels in block
        for ll, ilev in enumerate(en[1][ind]['ILEV']):
            # Dictionary to store data
            data[ilev] = {}

            # Level energy
            data[ilev]['E_eV'] = en[1][ind]['ENERGY'][ll]

            # Level qunatum numbers
            data[ilev]['VNL'] = en[1][ind]['VNL'][ll]

            # Level name
            data[ilev]['sname'] = en[1][ind]['sname'][ll]
            data[ilev]['name'] = en[1][ind]['name'][ll]

            # Number of electrons
            data[ilev]['nele'] = en[1][ind]['NELE']

            # If solving for this level population
            if ((int(data[ilev]['VNL']/100) <= nsim)
                #and (data[ilev]['E_eV'] > 0)
                and (data[ilev]['nele'] == nele)):
                lvls.append(ilev)

            # Formats data
            if mode == 'rt':
                data = _format_rt(
                    data = data,
                    ilev = ilev,
                    rt1 = rt1,
                    rt2 = rt2,
                    rt3 = rt3,
                    rt4 = rt4,
                    rt6 = rt6,
                    )
            elif mode == 'mr':
                data = _format_mr(
                    data = data,
                    ilev = ilev,
                    tr = tr,
                    ce = ce,
                    rr = rr,
                    ci = ci,
                    Te_eV = Te_eV,
                    )
       
    # Ouput
    return data, lvls

# Format FAC data generated with pfac.fac.MaxwellRate
def _format_mr(
    data = None,
    ilev = None,
    tr = None,
    ce = None,
    rr = None,
    ci = None,
    Te_eV = None,
    ):
    '''
    TO DO:
        1) 2-photon emission?
    '''

    # --- Finds the spontaneous emission rates for this state --- #
    data[ilev]['spon_emis'] = {} # [1/s]

    # Loop over blocks
    for blk in np.arange(len(tr[1])):
        # Loop over transitions in block
        for tran in np.arange(len(tr[1][blk]['upper_index'])):
            # If this is an upper state of the data
            if tr[1][blk]['upper_index'][tran] == ilev:
                # Stores data labeled by lower level
                lwr = tr[1][blk]['lower_index'][tran]
                data[ilev]['spon_emis'][lwr] = tr[1][blk]['rate'][tran]

    # --- Finds the collisional (de-)excitation for this state --- #
    data[ilev]['coll_excit'] = {} # [cm3/s]

    # Loop over upper levels
    if ilev in ce.keys():
        for upr in ce[ilev]['coll_excit'].keys():
            data[ilev]['coll_excit'][upr] = float(interp1d(
                ce['Te_eV'],
                ce[ilev]['coll_excit'][upr]
                )(Te_eV))

    # --- Finds the collisional dexcitation for this state --- #
    data[ilev]['coll_de_excit'] = {} # [cm3/s]

    # Loop over levels
    for lwr in ce.keys():
        if lwr != 'Te_eV':
            for upr in ce[lwr]['coll_de_excit'].keys():
                if upr == ilev:
                    data[ilev]['coll_de_excit'][lwr] = float(interp1d(
                        ce['Te_eV'],
                        ce[lwr]['coll_de_excit'][upr]
                        )(Te_eV))

    # --- Finds the radiative recombination for this state --- #
    data[ilev]['rad_recomb'] = {} # [cm3/s]

    # If this is a bound state
    if ilev in rr.keys():
        # Loop over continuum states
        for cnt in rr[ilev]['rad_recomb'].keys():
            data[ilev]['rad_recomb'][cnt] = float(interp1d(
                rr['Te_eV'],
                rr[ilev]['rad_recomb'][cnt]
                )(Te_eV))

    # If this is a continuum state
    else:
        # Loop over bound states
        for bnd in rr.keys():
            if bnd == 'Te_eV':
                continue
            if ilev in rr[bnd]['rad_recomb'].keys():
                # Also stores for continuum state
                data[ilev]['rad_recomb'][bnd] = float(interp1d(
                    rr['Te_eV'],
                    rr[bnd]['rad_recomb'][ilev]
                    )(Te_eV))

    # --- Finds the collisional ionization for this state --- #
    data[ilev]['coll_ion'] = {} # [cm3/s]

    # Loop over continuum levels
    if ilev in ci.keys():
        for cnt in ci[ilev]['coll_ion'].keys():
            data[ilev]['coll_ion'][cnt] = float(interp1d(
                ci['Te_eV'],
                ci[ilev]['coll_ion'][cnt]
                )(Te_eV))

    # Output
    return data

# Format FAC data generated with pfac.crm
def _format_rt(
    data = None,
    ilev = None,
    rt1 = None,
    rt2 = None,
    rt3 = None,
    rt4 = None,
    rt6 = None,
    ):

    # --- Finds the spontaneous emission rates for this state --- #
    data[ilev]['spon_emis'] = {} # [1/s]

    # Loop over transitions in rt1
    for tran in np.arange(rt1.shape[0]):
        # If this is an upper state of the data
        if rt1[tran,1] == ilev:
            # stores data labeled by lower level
            lwr = rt1[tran,2]
            data[ilev]['spon_emis'][lwr] = rt1[tran,3]

    # --- Finds the 2-photon spontaneous emission rates for this state --- #

    # Loop over transitions in rt2
    for tran in np.arange(rt2.shape[0]):
        # If this is an upper state of the data
        if rt2[tran,1] == ilev:
            # stores data labeled by lower level
            lwr = rt2[tran,2]
            data[ilev]['spon_emis'][lwr] += rt2[tran,3]

    # --- Finds the collisional excitation for this state --- #
    data[ilev]['coll_excit'] = {} # [cm3/s]

    # Loop over transitions in rt3
    for tran in np.arange(rt3.shape[0]):
        # If this is a lower level of the data
        if rt3[tran,1] == ilev:
            # stores data labeled by upper level
            upr = rt3[tran,2]
            data[ilev]['coll_excit'][upr] = rt3[tran,3]*1e-10

    # --- Finds the collisional de-excitation for this state --- #
    data[ilev]['coll_de_excit'] = {} # [cm3/s]

    # Loop over transitions in rt3
    for tran in np.arange(rt3.shape[0]):
        # If this is a upper level of the data
        if rt3[tran,2] == ilev:
            # stores data labeled by upper level
            lwr = rt3[tran,1]
            data[ilev]['coll_de_excit'][lwr] = rt3[tran,4]*1e-10

    # --- Finds the radiative recombination for this state --- #
    data[ilev]['rad_recomb'] = {} # [cm3/s]

    # Loop over transitions in rt3
    for tran in np.arange(rt4.shape[0]):
        # If this is a bound level of the data
        if rt4[tran,2] == ilev:
            # stores data labeled by coninuum level
            cnt = rt4[tran,1]
            data[ilev]['rad_recomb'][cnt] = rt4[tran,3]*1e-10

        # If this is a continuum level of the data
        elif rt4[tran,1] == ilev:
            # stores data labeled by bound level
            bnd = rt4[tran,2]
            data[ilev]['rad_recomb'][bnd] = rt4[tran,3]*1e-10

    # --- Finds the collisional ionization for this state --- #
    data[ilev]['coll_ion'] = {} # [cm3/s]

    # Loop over transitions in rt3
    for tran in np.arange(rt6.shape[0]):
        # If this is a bound level of the data
        if rt6[tran,1] == ilev:
            # stores data labeled by coninuum level
            cnt = rt6[tran,2]
            data[ilev]['coll_ion'][cnt] = rt6[tran,3]*1e-10

    # Output
    return data
 
# Reads necessary FAC data
def _read_FAC_mr(
    path_fac = None,
    ion = None,
    nele = None,
    ):

    # Loads energy level data
    en = rfac.read_en(
        os.path.join(
            path_fac, 
            ion+'%02d'%(nele)+'a.en'
            )
        )

    # Loads spontaneous emission data
    tr = rfac.read_tr(
        os.path.join(
            path_fac, 
            ion+'%02d'%(nele)+'a.tr'
            )
        )

    # Loads Maxwellian-averaged collisional (de-)excitation
    ce = _read_mr(
        path_fac  = path_fac,
        ion = ion,
        nele = nele,
        data = 'ce',
        )

    # Loads Maxwellian-averaged radiative recombination
    rr = _read_mr(
        path_fac  = path_fac,
        ion = ion,
        nele = nele,
        data = 'rr',
        )

    # Loads Maxwellian-averaged collisional ionization
    ci = _read_mr(
        path_fac  = path_fac,
        ion = ion,
        nele = nele,
        data = 'ci',
        )

    # Output
    return en, tr, ce, rr, ci

# Read Maxwellian-averaged data files
def _read_mr(
    path_fac = None,
    ion = None,
    nele = None,
    data = None,
    ):

    # Reads data file
    f = open(
        os.path.join(
            path_fac,
            ion+'%02d'%(nele)+data+'.mr'
            ),
        'r'
        )

    # Initializes output dictionary
    out = {}

    # Data orginaization labels
    if data == 'ce':
        label = 'coll_excit'
        incl_inv = True
        label_inv = 'coll_de_excit'
    elif data == 'rr':
        label = 'rad_recomb'
        incl_inv = False
    elif data == 'ci':
        label = 'coll_ion'
        incl_inv = False

    # Loop over lines
    for line in f:
        # Skip line breaks
        if line == '\n':
            continue

        # If reading a header
        if line.split(' ')[0] == '#':
            # Lower level
            lwr = int(line.split('\t')[0][1:])
            # Upper level
            upr = int(line.split('\t')[2])

            # Statistical weights
            lwr_2J = float(line.split('\t')[1])
            upr_2J = float(line.split('\t')[3])

            # Energy difference
            dE = float(line.split('\t')[4]) # [eV]

            # Number of temperature points
            ntemp = int(line.split('\t')[5])
            indt = 0

            # If temperature mesh not yet included
            if 'Te_eV' not in out.keys():
                out['Te_eV'] = []

            # If new lower state
            if lwr not in out.keys():
                out[lwr] = {}
                out[lwr][label] = {}
                if incl_inv:
                    out[lwr][label_inv] = {}

            out[lwr][label][upr] = []
            if incl_inv:
                out[lwr][label_inv][upr] = []

        # If reading data
        else:
            line = line.replace('\t', ' ')
            # If need to add temperature mesh
            if len(out['Te_eV']) < ntemp:
                out['Te_eV'].append(
                    float(line.split('  ')[0])
                    )

            # Adds rate coefficient data, [cm3/s]
            out[lwr][label][upr].append(
                float(line.split('  ')[-1])*1e-10
                )

            # If including the inverse process
            if incl_inv:
                out[lwr][label_inv][upr].append(
                    out[lwr][label][upr][-1]
                    *(lwr_2J+1)/(upr_2J+1)
                    *np.exp(dE/out['Te_eV'][indt])
                    )

                # Increases index to next temp point
                indt += 1

    # Output
    return out

# Reads necessary FAC data
def _read_FAC_rt(
    path_fac = None,
    ion = None,
    nele = None,
    ):

    # Loads energy level data
    en = rfac.read_en(
        os.path.join(
            path_fac, 
            ion+'%02d'%(nele)+'a.en'
            )
        )

    # Loads in rate data
    rt1 = np.loadtxt(
        os.path.join(
            path_fac, 
            ion+'%02d'%(nele)+'a.rt1'
            )
        ) # Spontaneous emissiion, dim(ntrans, 5), [1/s]

    rt2 = np.loadtxt(
        os.path.join(
            path_fac, 
            ion+'%02d'%(nele)+'a.rt2'
            )
        ) # 2-photon spontaneous emissiion, dim(ntrans, 5), [1/s]
    if rt2.ndim == 1:
        rt2 = rt2[None,:]

    rt3 = np.loadtxt(
        os.path.join(
            path_fac, 
            ion+'%02d'%(nele)+'a.rt3'
            )
        ) # Collisional (de)excitation, dim(ntrans, 5), [1e-10 cm3/s]

    rt4 = np.loadtxt(
        os.path.join(
            path_fac, 
            ion+'%02d'%(nele)+'a.rt4'
            )
        ) # Radiative recombination, dim(ntrans, 5), [1e-10 cm3/s]

    rt6 = np.loadtxt(
        os.path.join(
            path_fac, 
            ion+'%02d'%(nele)+'a.rt6'
            )
        ) # Collisional ionization, dim(ntrans, 5), [1e-10 cm3/s]

    # Output
    return en, rt1, rt2, rt3, rt4, rt6


#########################################
#
#           Terms
#
#########################################

'''

Population equation into state rr of charge Z

Assume:
    1) Optically thin 
        -> Ignore radiation field interactions
    2) Low energy density 
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
    + src_ion                               # Source term,
                                                # ionizing charge Z-1
    + src_rec                               # Source term,
                                                # recombining charge Z+1
    
sink = -1*(
    \sum_{E_ll < E_rr} A_{rr->ll}       # Decay to lower states
    + N_e * S_{rr -> cnt}               # Collisional ionization to continuum
    + N_e * \sum{ll=/=rr} S_{rr->ll}    # Collisional excitation from state rr
    )                                   # [1/s]

src_lwr = N_e * S_{ll->rr}^             # Collisional excitation (excl. ground) to state rr
                                        # Incl. Cascasde from unsimulated E_kk > E_rr
                                        # [1/s]

src_upr = (
    N_e * S_{ll->rr}                    # Collisional excitation to state rr
    + A_{ll->rr}                        # Spontaneous emission to state rr
    )                                   # [1/s]

src_ion = N_e * N_gg^{Z-1} * S_{gg->rr} # Collisional ionization from Z-1 ground
                                        # [1/cm3/s]

src_rec = (
    N_e * N_gg^{Z+1} * ar_{gg->rr}^     # Radiative recombination from Z+1 ground
                                        # Incl. Cascade from unsimulated E_kk > E_rr
    + N_e * N_gg^{Z+1} * ad_{gg->rr}^   # Dielectronic capture from Z+1 ground
                                        # Incl. Cascade from unsimulated E_kk > E_rr
    )                                   # [1/cm3/s]

'''

# Finds cascades from ignorably popullated states into state upr
def _calc_cascade(
    quant = 'coll_excit',
    data_all = None,
    # Levels in consideration
    lwr = None,
    upr = None,
    # Plasma
    ne_cm3 = None,
    nsim = None,
    ):

    # Obtains rate coefficient data from lower state, [cm3/s]
    data = data_all[lwr][quant]

    # Initializes cascade term
    casc = 0.0

    # Loop over lvls
    for higher in data.keys():
        # If higher energy state
        if ((data_all[upr]['E_eV'] < data_all[higher]['E_eV'])
            and (int(data_all[higher]['VNL']/100) > nsim)):
            # Gets data from lower to higher
            data_hu = data[higher]

            # Calculates branching ratio from higher to upper
            try:
                b_hu_top = (
                    data_all[higher]['spon_emis'][upr]
                    + ne_cm3 * data_all[higher]['coll_de_excit'][upr]
                    )
                b_hu_bot  = 0
                for key in data_all[higher]['spon_emis'].keys():
                    b_hu_bot += data_all[higher]['spon_emis'][key]
                for key in data_all[higher]['coll_excit'].keys():
                    b_hu_bot += ne_cm3 *data_all[higher]['coll_excit'][key]
                for key in data_all[higher]['coll_de_excit'].keys():
                    b_hu_bot += ne_cm3 *data_all[higher]['coll_de_excit'][key]
                for key in data_all[higher]['coll_ion'].keys():
                    b_hu_bot += ne_cm3 *data_all[higher]['coll_ion'][key]

                casc += data_hu *b_hu_top/b_hu_bot
                
            # In case there's no coupling 
            except:
                print('skip cascade')
                print('quantity= '+quant)
                print('lower= '+ str(int(lwr)))
                print('higher= '+ str(int(higher)))
                print('upper='+ str(int(upr)))
                print('\n')

    # [1/s]
    return casc *ne_cm3

# Calculates the sink term in population equation
def _calc_sink(
    data = None,
    state = None,
    ne_cm3 = None,
    ):

    # Init
    sink = 0.0

    # Radiative decay out of state rr
    for other in data[state]['spon_emis'].keys():
        sink -= data[state]['spon_emis'][other]

    # Collisional excitation/de-excitation out of state rr
    for other in data[state]['coll_excit'].keys():
        sink -= ne_cm3 *data[state]['coll_excit'][other]
    for other in data[state]['coll_de_excit'].keys():
        sink -= ne_cm3 *data[state]['coll_de_excit'][other]

    # Collisional ionization out of state rr
    for other in data[state]['coll_ion'].keys():
        sink -= ne_cm3 *data[state]['coll_ion'][other]

    # [1/s]
    return sink

# Source from lower energy states
def _calc_src_lwr(
    data = None,
    lwr = None,
    upr = None,
    ne_cm3 = None,
    nsim = None,
    ):

    # Init
    src = 0.0

    # Collisional excitation from cc into rr including cascades
    try:
        src += ne_cm3 *data[lwr]['coll_excit'][upr]
    except:
        src += 0.0

    src += _calc_cascade(
        quant = 'coll_excit',
        data_all = data,
        lwr = lwr,
        upr = upr,
        ne_cm3=ne_cm3,
        nsim = nsim,
        )

    # [1/s]
    return src

# Source from higher energy states
def _calc_src_upr(
    data = None,
    upr = None,
    lwr = None,
    ne_cm3 = None,
    ):

    # Init
    src = 0.0

    # Collisional de-excitation from cc to rr
    try:
        src += ne_cm3 *data[upr]['coll_de_excit'][lwr]
    except:
        src += 0.0

    # Radiative decay from cc to rr
    try:
        src += data[upr]['spon_emis'][lwr]
    except:
        src += 0.0

    # [1/s]
    return src


# Source from charge Z+1 ground state recombining to state rr
def _calc_src_rec(
    data = None,
    state = None,
    ne_cm3 = None,
    nz_cm3 = None,
    ZZ = None,
    nsim = None,
    ):
    '''
    TO DO:
        1) Add dielectronic recombination
    '''

    # Init 
    src = 0.0

    # Ground state of ion Z+1 radiative recombines to rr including cascades
    for rec in data[state]['rad_recomb'].keys():
        src += ne_cm3 *data[state]['rad_recomb'][rec] * nz_cm3[ZZ +1]

        src += _calc_cascade(
            quant = 'rad_recomb',
            data_all = data,
            lwr = rec,
            upr = state,
            ne_cm3=ne_cm3,
            nsim=nsim,
            ) *nz_cm3[ZZ +1]

    # [1/cm3/s]
    return src

# Source from charge Z+1 ground state ionizing to state rr
def _calc_src_ion(
    data = None,
    ):
    '''
    TO DO:
        1) Add collisional ionization
    '''

    # Init 
    src = 0.0

    # [1/cm3/s]
    return src


