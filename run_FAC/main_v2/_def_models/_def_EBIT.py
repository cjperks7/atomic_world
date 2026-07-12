'''

Stores what levels were used to model the EBIT spectra

cjperks
July 12th, 2026

'''

import numpy as np

__all__ = [
    '_EBIT_W_Lshell',
    ]

########################################################
#
#       Ground states
#
########################################################

# W L-shell from Ne- to Ar-like
def _EBIT_W_Lshell(
    nouter = 0
    ):
    llist = [
        's', 'p', 'd', 'f', 
        'g', 'h', 'i',  # skip 'j'
        'k', 'l', 'm'
        ]   # max at ll=(nn-1); 0, 1, 2, 3,     4, 5, 6,     7, 8, 9,

    # Number of levels in ground state
    ns_orb = min(2, nouter)
    np_orb = min(6, max(nouter-2,0))

    dout = {}

    ########################################################
    #
    #       Settings
    #
    ########################################################

    # Define energy levels by a direct string
    dout['method'] = 'string' 

    dout['spectator'] = {           # --- Highly excited (spectator) settings --- #
        'ion_ind': 0,       # Ionized state to base spectators off of
        'rec_pw_max': 3,    # Maximum orbital angular momentum
        'n_frozen': 7,      # n>n_frozen uses frozen core approximation
        'n_spec_min': 8,    # Minimum quantum number of spectator
        'n_spec_max': 7,    # Maximum quantum number of spectator
        'n_ext': None,      # Highly excited state if desired
        }
    dout['restrictions'] = {
        'skip_CI': True,    # Skip configuration interactions (ce: exc->exc)
        }

    ########################################################
    #
    #       Ground states
    #
    ########################################################
    dout['grd'] = {}

    kkey = 'grd.0'

    if nouter == 0:
        ktrans = '1*2 2*8'
    elif nouter <= 2:
        ktrans = '1*2 2*8 3*%i;3p<2;3d<2'%(nouter)
    else:
        ktrans = '1*2 2*8 3*%i;3p<%i;3d<2'%(nouter, min(nouter,7))

    # Ignore exciting both 3s electrons
    if nouter >1:
        ktrans += ';3s>0'

    dout['grd'][kkey] = ktrans

    ########################################################
    #
    #       Excited states
    #
    ########################################################
    dout['exc'] = {}

    ### --- 3l --> nl', outershell transitions --- ###
    dlims = {
        0: {'nmin': 4, 'nmax': 4, 'lmin': 0, 'lmax': 3},    # Ne-like, irrelavent here
        1: {'nmin': 4, 'nmax': 4, 'lmin': 0, 'lmax': 3},    # Na-like
        2: {'nmin': 4, 'nmax': 4, 'lmin': 0, 'lmax': 3},    # Mg-like
        3: {'nmin': 4, 'nmax': 4, 'lmin': 0, 'lmax': 3},    # Al-like
        4: {'nmin': 4, 'nmax': 4, 'lmin': 0, 'lmax': 3},    # Si-like
        5: {'nmin': 4, 'nmax': 4, 'lmin': 0, 'lmax': 3},    # P-like
        6: {'nmin': 4, 'nmax': 4, 'lmin': 0, 'lmax': 3},    # S-like
        7: {'nmin': 4, 'nmax': 4, 'lmin': 0, 'lmax': 3},    # Cl-like
        8: {'nmin': 4, 'nmax': 4, 'lmin': 0, 'lmax': 3},    # Ar-like
        }
    nmin = dlims[nouter]['nmin']    # NOTE: n=3 levels handled by ground state
    nmax = dlims[nouter]['nmax']    
    lmin = dlims[nouter]['lmin']    
    lmax = dlims[nouter]['lmax']    

    kshells = ['3s', '3p']

    # Loop over outershells
    for kshell in kshells:
        # If exciting 3s electron, but no 3s electron
        if kshell == '3s' and nouter <=0:
            continue
        elif kshell == '3s':
            ns_exc = 1
            np_exc = 0
        # If exciting 3p electron, but no 3p electron
        if kshell == '3p' and nouter <=2:
            continue
        elif kshell == '3p':
            ns_exc = 0
            np_exc = 1

        # Loop over n
        for nn in np.arange(nmin, nmax+1):
            # Loop over l
            for ll in np.arange(lmin, lmax+1):
                # Skip nonsense
                if ll >= nn:
                    continue

                # Labels transition
                kkey = 'exc.%s->%i%s'%(
                    kshell, nn, llist[ll]
                    )

                # Init
                ktrans = '1*2 2*8'

                # Handles 3s-> nl, n>3
                ktrans += (
                    (
                        ' 3s%i'%(ns_orb-ns_exc) if ns_orb>ns_exc
                        else ''
                        )
                    + (
                        ' 3p%i'%(np_orb-np_exc) if np_orb>np_exc
                        else ''
                        )
                    + ' %i%s1'%(nn, llist[ll])
                    )

                # Makes dictionary entry
                dout['exc'][kkey] = ktrans

    ### --- 2p -> nl, innershell transitions --- ###
    dlims = {
        0: {'nmin': 3, 'nmax': 4, 'lmin': 0, 'lmax': 3},    # Ne-like
        1: {'nmin': 3, 'nmax': 4, 'lmin': 0, 'lmax': 3},    # Na-like
        2: {'nmin': 3, 'nmax': 4, 'lmin': 0, 'lmax': 3},    # Mg-like
        3: {'nmin': 3, 'nmax': 4, 'lmin': 0, 'lmax': 3},    # Al-like
        4: {'nmin': 3, 'nmax': 4, 'lmin': 0, 'lmax': 3},    # Si-like
        5: {'nmin': 3, 'nmax': 4, 'lmin': 0, 'lmax': 3},    # P-like
        6: {'nmin': 3, 'nmax': 4, 'lmin': 0, 'lmax': 3},    # S-like
        7: {'nmin': 3, 'nmax': 4, 'lmin': 0, 'lmax': 3},    # Cl-like
        8: {'nmin': 3, 'nmax': 4, 'lmin': 0, 'lmax': 3},    # Ar-like
        }
    nmin = dlims[nouter]['nmin']  
    nmax = dlims[nouter]['nmax']    
    lmin = dlims[nouter]['lmin']    
    lmax = dlims[nouter]['lmax']  

    # Loop over n
    for nn in np.arange(nmin, nmax+1):
        # Loop over l
        for ll in np.arange(lmin, lmax+1):
            # If 3s orbital is filled
            if nn==3 and ll==0 and nouter>=2:
                continue
            # If 3p orbital is filled
            if nn==3 and ll==1 and nouter>=8:
                continue
            if ll >= nn:
                continue
            
            # Labels transition
            kkey = 'exc.2p->%i%s'%(
                nn, llist[ll]
                )

            # Init
            ktrans = '1*2 2s2 2p5'

            # Handles 2p->3s orbital
            if nn==3 and ll==0:
                ktrans += (
                    ' 3s%i'%(nouter+1)
                    )

            # Handles 2p->3p orbital
            elif nn==3 and ll==1:
                ktrans += (
                    (
                        ' 3s%i'%(min(2,nouter)) if nouter >0
                        else ''
                        )
                    + (
                        ' 3p1' if nouter < 3
                        else ' 3p%i'%(nouter-1)
                        )
                    )

            # Handles 2p-> nl, n>3
            else:
                ktrans += (
                    (
                        ' 3s%i'%(ns_orb) if ns_orb>0
                        else ''
                        )
                    + (
                        ' 3p%i'%(np_orb) if np_orb>0
                        else ''
                        )
                    + ' %i%s1'%(nn, llist[ll])
                    )

            # Makes dictionary entry
            dout['exc'][kkey] = ktrans

    ########################################################
    #
    #       Ground states
    #
    ########################################################
    dout['ion'] = {}

    ### --- Outshell ionization --- ###
    kkey = 'ion.0'

    if nouter ==0:
        ktrans = '1*2 2s2 2p5'
    elif nouter==1:
        ktrans = '1*2 2*8'
    elif nouter <=3:
        ktrans = '1*2 2*8 3s%i'%(nouter-1)
    else:
        ktrans = '1*2 2*8 3s2 3p%i'%(nouter-3)

    # Makes dictionary entry
    dout['ion'][kkey] = ktrans

    ### --- Innershell ionization --- ###
    if nouter >0:
        kkey = 'ion.1'

        if nouter <=2:
            ktrans = '1*2 2s2 2p5 3s%i'%(nouter)
        else:
            ktrans = '1*2 2s2 2p5 3s2 3p%i'%(nouter-2)

        dout['ion'][kkey] = ktrans

    # Output
    return dout
