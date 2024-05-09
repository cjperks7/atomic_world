'''

_ce_tr.py is a script meant to facilitate calculating
collisional ionization and radiative recombination using FAC

cjperks
May 9, 2024

'''

# Modules
import numpy as np
import sys, os
import copy
import time

__all__ = [
    '_run_rr',
    '_run_ci'
    ]


######################################################
#
#              Radiative recombination
#
######################################################

# Radiative recombination rates
def _run_rr(
    name = None,
    groups = None,
    settings = None,
    fac = None,
    ):

    # Track time spent
    start = time.time()
    start0 = start

    print('RR: Calc radiative recombination')

    # Loop over types
    for typ1 in groups.keys():
        # Loop over states
        for ind1, st1 in enumerate(groups[typ1]):

            # Loop over ionized states
            for st2 in groups['ion']:

                # Recombining to spectator state
                if typ1 == 'ion' and ind1 == settings['spectator']['ion_ind']:
                    st1 = groups[typ1][settings['spectator']['ion_ind']]

                    # Loop over spectator electrons
                    for nn1 in range(
                        settings['spectator']['n_spec_min'], 
                        settings['spectator']['n_spec_max']+1
                        ):
                        _calc_rr(
                            fac=fac, name=name, ind=ind1,
                            st1 = ([st1], nn1),
                            st2 = [st2]
                            )
                    # Highly excited state
                    if settings['spectator']['n_ext'] is not None:
                        _calc_rr(
                            fac=fac, name=name, ind=ind1,
                            st1 = ([st1], settings['spectator']['n_ext']),
                            st2 = [st2]
                            )

                # Recombining to ground or excitated state
                elif typ1 in ['grd','exc']:
                    _calc_rr(fac=fac, name=name, ind=ind1, st1=[st1], st2=[st2])

    # Print the time this took
    stop = time.time()
    s='Done %10.3E s'%(stop-start)
    print(s)

# Calculates radiative recombination
def _calc_rr(
    fac = None,
    name = None,
    ind = None,
    st1 = None,
    st2 = None,
    ):

    # Calculates rad. recomb and photoionization rates
    print(str(st1)+' --> '+str(st2))
    fac.RRTable(name+'b.rr', st1, st2)
    fac.SetPEGrid(0)
    fac.SetUsrPEGrid(0)
    fac.SetRRTEGrid(0)

    # Clears memory
    if ind > 0:
        fac.Reinit(radial = 1, recombination = 1, ionization = 1)

######################################################
#
#              Collisional ionization
#
######################################################

# Collisional ionization cross-section
def _run_ci(
    name = None,
    groups = None,
    settings = None,
    fac = None,
    ):

    # Track time spent
    start = time.time()
    start0 = start

    print('CI: Calc collisional ionization')

    # Loop over types
    for typ1 in groups.keys():
        # Loop over states
        for ind1, st1 in enumerate(groups[typ1]):

            # Loop over ionized states
            for st2 in groups['ion']:

                # Ionizing from spectator state
                if typ1 == 'ion' and ind1 == settings['spectator']['ion_ind']:
                    st1 = groups[typ1][settings['spectator']['ion_ind']]

                    # Loop over spectator electrons
                    for nn1 in range(
                        settings['spectator']['n_spec_min'], 
                        settings['spectator']['n_spec_max']+1
                        ):
                        _calc_ci(
                            fac=fac, name=name, ind=ind1,
                            st1 = ([st1], nn1),
                            st2 = [st2]
                            )
                    # Highly excited state
                    if settings['spectator']['n_ext'] is not None:
                        _calc_ci(
                            fac=fac, name=name, ind=ind1,
                            st1 = ([st1], settings['spectator']['n_ext']),
                            st2 = [st2]
                            )

                # Ionizing from ground or excitated state
                elif typ1 in ['grd','exc']:
                    _calc_ci(fac=fac, name=name, ind=ind1, st1=[st1], st2=[st2])

    # Print the time this took
    stop = time.time()
    s='Done %10.3E s'%(stop-start)
    print(s)

# Calculates radiative recombination
def _calc_ci(
    fac = None,
    name = None,
    ind = None,
    st1 = None,
    st2 = None,
    ):

    # Calculates collisional excitation cross-sections
    print(str(st1)+' --> '+str(st2))
    fac.CITable(name+'b.ci', st1, st2)

    # Clears memory
    if ind > 0:
        fac.Reinit(radial = 1, recombination = 1, ionization = 1)