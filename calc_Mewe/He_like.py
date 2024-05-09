'''

Script to calculate He-like Kalpha PECs using the formulas in
R. Mewe and J. Schrijver, Astronomy & Astrophysics 65, 99-114 (1978)

cjperks
Jan 19th, 2024

'''

# Modules
import numpy as np


##################################################
#
#            Main
#
##################################################

def calc_He(
    sp = None,  # species
    Te = None,  # [eV], dim(ntemp,)
    incl_physics = 'all',
    ):

    # Shape of dim=0 is     [w,  1', x, y, 0, z]
    # In Mewe's notation    [m', 1', 2, 1, 0, m]

    if incl_physics == 'all':
        incl_physics = ['ex', 'rr', 'dr']

    # init output
    out = {}

    # Loop over physics
    for pp in incl_physics:
        # Collisional excitation
        if pp == 'ex':
            out['ex'] = {}
            (
                out['ex']['PEC'], 
                out['ex']['cas'],
                out['ex']['Omega'],
                out['ex']['Egrid']
                ) = _get_ex(
                sp=sp,
                Te=Te,
                )

        # Radiative recombination
        elif pp == 'rr':
            out['rr'] = {}
            out['rr']['PEC'], out['rr']['cas'] = _get_rr(
                sp=sp,
                Te=Te,
                )

        # Dielectronic recombination
        elif pp == 'dr':
            out['dr'] = {}
            out['dr']['PEC'] = _get_dr(
                sp=sp,
                Te=Te,
                )

    # Output
    return out

##################################################
#
#           Utilities
#
##################################################

# Calculates dielectronic recombination coefficient, including cascades
def _get_dr(
    sp = None,
    Te = None, # [eV], dim(ntemp,)
    ):

    # Init output
    ad = np.zeros((6, len(Te)))

    # Nuclear charge
    if sp == 'Ar':
        ZZ = 18
        wk = [5,3,1]

    # Eqn 62
    pp = 33 *(ZZ-1)**0.6 *Te**-0.3

    # z-line, Eqn 63
    ad[-1,:] = 6.46e-8 *ZZ**4 *Te**(-3/2) *(
        9/(1 +7e-5 *ZZ**4) *np.exp(-78900 *(ZZ+0.5)**2/Te)
        + 27/(1 +8e-5 *ZZ**4) *np.exp(-101800 *ZZ**2/Te)
        + 380 *(1+pp)**-1/(1 +5e-3 *ZZ**3) *np.exp(-118400 *ZZ**2/Te)
        )

    # w-line, Eqn 64
    ad[0,:] = 6.46e-8 *ZZ**4 *Te**(-3/2) *(
        3/(1 +3e-6 *ZZ**4) *np.exp(-78900 *(ZZ+0.5)**2/Te)
        + 0.5/(1 +2.2e-5 *ZZ**4) *np.exp(-101800 *ZZ**2/Te)
        + 6.3/(1 +5e-3 *ZZ**3) *np.exp(-118400 *ZZ**2/Te)
        )

    # x-line, Eqn 65
    ad[2,:] = 6.46e-8 *(wk[0]/9) *ZZ**4 *Te**(-3/2) *(
        18/9.5 *np.exp(-78900 *(ZZ+0.5)**2/Te)
        + 54/(1 +1.9e-4 *ZZ**4) *np.exp(-101800 *ZZ**2/Te)
        + 380/(1 +5e-3 *ZZ**3) *pp/(1+pp) *np.exp(-118400 *ZZ**2/Te)
        )

    # y-line, Eqn 65
    ad[3,:] = 6.46e-8 *(wk[1]/9) *ZZ**4 *Te**(-3/2) *(
        18/9.5 *np.exp(-78900 *(ZZ+0.5)**2/Te)
        + 54/(1 +1.9e-4 *ZZ**4) *np.exp(-101800 *ZZ**2/Te)
        + 380/(1 +5e-3 *ZZ**3) *pp/(1+pp) *np.exp(-118400 *ZZ**2/Te)
        )

    # 0-line, Eqn 65
    ad[4,:] = 6.46e-8 *(wk[2]/9) *ZZ**4 *Te**(-3/2) *(
        18/9.5 *np.exp(-78900 *(ZZ+0.5)**2/Te)
        + 54/(1 +1.9e-4 *ZZ**4) *np.exp(-101800 *ZZ**2/Te)
        + 380/(1 +5e-3 *ZZ**3) *pp/(1+pp) *np.exp(-118400 *ZZ**2/Te)
        )

    # 1'-line
    ad[1,:] = 6.46e-8 *ZZ**4 *Te**(-3/2) *(
        12/(1 +6e-6 *ZZ**4) *np.exp(-78900 *(ZZ+0.5)**2/Te)
        + 18/(1 +3e-5 *ZZ**4) *np.exp(-101800 *ZZ**2/Te)
        + 69/(1 +5e-3 *ZZ**3) *np.exp(-118400 *ZZ**2/Te)
        )

    # Output, [cm3/s], dim(6, ntemp)
    return ad

# Calculates radiative recombination coefficient, including cascades
def _get_rr(
    sp = None,
    Te = None, # [eV], dim(ntemp,)
    ):

    # Init output
    ar = np.zeros((6, len(Te)))
    dar = np.zeros((6, len(Te)))

    # Charge
    if sp == 'Ar':
        zz = 18-2

    # z-line, Eqn 56
    ar[-1,:] = 5.9e-18 * zz**1.8 * Te**-0.4
    dar[-1,:] = 1e-11 *zz**2.2 *Te**-0.6

    # x-line, Eqn 57
    ar[2,:] = 3.6e-11 *(5/9) *zz**2.4 *Te**-0.7
    dar[2,:] = 3.6e-10 *(5/9) *zz**2.8 *Te**-0.9

    # y-line
    ar[3,:] = 3.6e-11 *(3/9) *zz**2.4 *Te**-0.7
    dar[3,:] = 3.6e-10 *(3/9) *zz**2.8 *Te**-0.9

    # 0-line
    ar[4,:] = 3.6e-11 *(1/9) *zz**2.4 *Te**-0.7
    dar[4,:] = 3.6e-10 *(1/9) *zz**2.8 *Te**-0.9

    # w-line, Eqn 58
    ar[0,:] = 5.4e-13 *zz**2 *Te**-0.5
    dar[0,:] = ar[0,:] *dar[-1,:]/ar[-1,:]

    # 1'-line, Eqn 59
    ar[1,:] = 1.4e-11 *zz**2.52 *Te**-0.76
    dar[1,:] = ar[1,:] *dar[2,:]/ar[2,:]

    # Output
    return ar+dar, 1+dar/ar

# Calcualtes excitation rate coefficient, including cascades
def _get_ex(
    sp = None,
    Te = None, # [eV], dim(ntemp,)
    ):

    # Calculate transition energy
    E_gk, Z_gk = _get_Egk(
        sp=sp,
        ) # [eV], dim(6,)
    yy = E_gk[:,None]/Te[None,:] # dim(6, ntemp)

    # Calculates Maxwell-averaged collision strength
    Omega_bar, Omega, Egrid = _get_omega(
        sp=sp,
        E_gk=E_gk,
        yy=yy,
        Z_gk=Z_gk
        ) # dim(6,ntemp)

    # Calculates excitation rate coefficient, Eqn 33
    S_gk = (
        8.62e-6 * Omega_bar
        * Te[None,:]**(-1/2)
        * np.exp(-yy)
        ) # [cm3/s], dim(6, ntemp)

    # Includes cascade effects
    S_gk_cas = S_gk.copy()

    # Get the ionization potential
    chi2 = _get_chi(
        sp=sp,
        ii = 2,
        )

    # For the z-line, Eqn 38
    dE = 0.16 *chi2
    S_gk_cas[-1,:] += 0.4*(
        S_gk[2,:]+S_gk[3,:]+S_gk[4,:]
        )*np.exp(-dE/Te)

    # The rest, Eqn 41
    S_gk_cas[0,:] *= 1.1
    S_gk_cas[1,:] *= 1.0
    S_gk_cas[2,:] *= 1.065
    S_gk_cas[3,:] *= 1.065
    S_gk_cas[4,:] *= 1.065

    # Output, [cm3/s], dim(6,ntemp)
    return S_gk_cas, S_gk_cas/S_gk, Omega, Egrid

# Formula for the ionization potential
def _get_chi(
    sp =None,
    ii = None, # number of electrons
    ):

    nz = np.asarray([1,1,2,2,2])
    xi = np.asarray([1,2,1,2,3])

    # Table 1
    if sp == 'Ar':
        aa = np.asarray([
            -0.04,
            0.593,
            1.568,
            2.139,
            2.908,
            0.801
            ])
        Z = 18

    # Output, Eqn 9
    return (
        13.6 * nz[ii-1]**-2 *(Z-aa[ii-1])**2
        )


# Formula for Maxwell-averaged collision strength
def _get_omega(
    sp = None,
    E_gk=None, # dim(6,)
    yy=None, # dim(6, ntemp)
    Z_gk=None, # dim(6,)
    ):

    if sp == 'Ar':
        wk = [5, 3, 1]

    # Table 4
    Z2A = np.asarray([
        0.8,    # m'
        0.6,    # 1'
        0.0,    # 2
        0.0,    # 1
        0.0,    # 0
        0.0     # m
        ])
    Z2B = np.asarray([
        -0.4,   # m'
        0.8,    # 1'
        0.0,    # 2
        0.0,    # 1
        0.0,    # 0
        0.0,    # m
        ])
    Z2C = np.asarray([
        0.0,            # m'
        0.0,            # 1'
        1.9*(wk[0]/9),  # 2
        1.9*(wk[1]/9),  # 1
        1.9*(wk[2]/9),  # 0
        0.6,            # m
        ])
    Z2D = np.asarray([
        0.0,    # m'
        0.0,    # 1'
        -0.25,  # 2
        -0.25,  # 1
        -0.25,  # 0
        -0.2,   # m
        ])
    Z2E = np.asarray([
        0.0,    # m'
        4.0,    # 1'
        0.0,    # 2
        0.0,    # 1
        0.0,    # 0
        0.0,    # m
        ])

    # Eqn 10
    aa = np.zeros(yy.shape) # dim(6, ntemp)
    aa[yy >= 1] = 0.5
    aa[yy < 1] = -0.5
    fy = (
        np.log((yy+1)/yy)
        - (
            0.36 + 0.03*(yy+0.01)**aa
            )
        * (yy+1)**-2
        ) # dim(6, ntemp)

    # Calculates collision strength, Eqn 35
    # NOTE: Incident electron energy
    Egrid = np.logspace(
        np.log10(E_gk),
        np.log10(E_gk*1e3),
        int(1e4),
        axis=1
        ) # [eV], dim(6,nE)
    UU = E_gk[:,None]/Egrid
    Omega = (
        Z2A[:,None]
        + Z2B[:,None] *UU**-1
        + Z2C[:,None] *UU**-2
        + 2*Z2D[:,None] *UU**-3
        + Z2E[:,None] *np.log(UU)
        )/Z_gk[:,None]**2

    # Calculates Maxwell-averaged collision strength, Eqn 36
    Omega_bar = (
        Z2A[:,None]
        + (
            Z2B[:,None]*yy
            - Z2C[:,None]*yy**2
            + Z2D[:,None]*yy**3
            + Z2E[:,None]
            )*fy
        + (Z2C[:,None] + Z2D[:,None])*yy
        - Z2D[:,None]*yy**2
        )/Z_gk[:,None]**2

    # Output
    return Omega_bar, Omega, Egrid

# Formula for transition energy
def _get_Egk(
    sp=None,
    ):

    # Table 3, screening factors
    if sp == 'Ar':
        a_gk = np.asarray([
            0.498,
            0.454,
            0.496,
            0.496,
            0.496,
            0.553
            ])
        Z = 18

    # Output, [eV], Eqn 37
    return 0.75 *13.6 *(Z -a_gk)**2, (Z-a_gk)