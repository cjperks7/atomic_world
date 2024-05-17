'''

_ce_tr.py is a script meant to facilitate calculating
collisional excitation and spontaneous emission using FAC

cjperks
May 9, 2024

'''

# Modules
import numpy as np
import sys, os
import copy
import time

__all__ = [
    '_run_tr',
    '_run_ce'
    ]


######################################################
#
#              Radiative transitions
#
######################################################

# Radiative transition rates
def _run_tr(
    name = None,
    groups = None,
    settings = None,
    fac = None,
    ):

    # Track time spent
    start = time.time()
    start0 = start

    print('TR: Calc radiative transtion')

    # List of state types
    types_tmp = list(groups.keys())

    # Loop over types
    for typ1 in groups.keys():
        # Loop over states in typ1
        for ind1, st1 in enumerate(groups[typ1]):

            # Loop over types
            for typ2 in types_tmp:
                # Loop over states in typ2
                for ind2, st2 in enumerate(groups[typ2]):

                    # If a ground state ...
                    if typ1 == 'grd':
                        # ... is interacting with a spectator state
                        if typ2 == 'ion' and ind2 == settings['spectator']['ion_ind']:
                            st2 = groups[typ2][settings['spectator']['ion_ind']]
                            for nn2 in range(
                                settings['spectator']['n_spec_min'], 
                                settings['spectator']['n_spec_max']+1
                                ):
                                _calc_tr(
                                    fac=fac, name = name,
                                    st1 = [st1],
                                    st2 = ([st2], nn2)
                                    )
                            # Highly excited state
                            if settings['spectator']['n_ext'] is not None:
                                _calc_tr(
                                    fac=fac, name=name,
                                    st1 = [st1],
                                    st2 = ([st2], settings['spectator']['n_ext'])
                                    )

                        # ... is interaction with a ground or excited state
                        elif typ2 in ['grd','exc']:
                            _calc_tr(fac=fac, name=name, st1=[st1], st2=[st2])

                    # If an excitation state ...
                    elif typ1 == 'exc':
                        # ... is interacting with an excited state
                        if typ2 == 'exc':
                            # Considers just higher or equal quantum number
                            if ind1 <= ind2:
                                _calc_tr(fac=fac, name=name, st1=[st1], st2=[st2])

                        # ... is interacting with a spectator state
                        elif typ2 == 'ion' and ind2 == settings['spectator']['ion_ind']:
                            st2 = groups[typ2][settings['spectator']['ion_ind']]
                            for nn2 in range(
                                settings['spectator']['n_spec_min'], 
                                settings['spectator']['n_spec_max']+1
                                ):
                                _calc_tr(
                                    fac=fac, name = name,
                                    st1 = [st1],
                                    st2 = ([st2], nn2)
                                    )
                            # Highly excited state
                            if settings['spectator']['n_ext'] is not None:
                                _calc_tr(
                                    fac=fac, name=name,
                                    st1 = [st1],
                                    st2 = ([st2], settings['spectator']['n_ext'])
                                    )
                        
                    # If a spectator state ...
                    elif typ1 == 'ion' and ind1 == settings['spectator']['ion_ind']:
                        st1 = groups[typ1][settings['spectator']['ion_ind']]
                        # ... is interacting with a spectator state
                        if typ2 == 'ion' and ind2 == settings['spectator']['ion_ind']:
                            st2 = groups[typ2][settings['spectator']['ion_ind']]
                            for nn1 in range(
                                settings['spectator']['n_spec_min'], 
                                settings['spectator']['n_spec_max']+1
                                ):
                                for nn2 in range(
                                    settings['spectator']['n_spec_min'], 
                                    settings['spectator']['n_spec_max']+1
                                    ):
                                    _calc_tr(
                                        fac=fac, name = name,
                                        st1 = ([st1], nn1),
                                        st2 = ([st2], nn2)
                                        )
                                # Highly excited state
                                if settings['spectator']['n_ext'] is not None:
                                    _calc_tr(
                                        fac=fac, name=name,
                                        st1 = ([st1], nn1),
                                        st2 = ([st2], settings['spectator']['n_ext'])
                                        )
                            # Highly excited state
                            if settings['spectator']['n_ext'] is not None:
                                _calc_tr(
                                    fac=fac, name=name,
                                    st1 = ([st1], settings['spectator']['n_ext']),
                                    st2 = ([st2], settings['spectator']['n_ext'])
                                    )

        # Makes sure we don't improperly order the energy levels
        types_tmp.remove(typ1)

    # Print the time this took
    stop = time.time()
    s='Done %10.3E s'%(stop-start)
    print(s)


# Calculates transition rates
def _calc_tr(
    fac = None,
    name = None,
    st1 = None,
    st2 = None,
    ):

    # Calculates radiative transition rates
    print(str(st1)+' --> '+str(st2) + ' 0')
    fac.TRTable(name+'b.tr',st1, st2, 0)
    # Clears memory
    fac.Reinit(radial = 1, excitation = 1)    


######################################################
#
#              Collisional excitation
#
######################################################


# Collisional excitation cross-sections
def _run_ce(
    name=None,
    groups=None,
    settings = None,
    fac = None
    ):

    # Limits for excited states to consider
    lim1 = 0
    lim2 = 0

    # Track time spent
    start = time.time()
    start0 = start

    print('CE: Calc collisional excitation')

    # List of state types
    types_tmp = list(groups.keys())

    # Loop over types
    for typ1 in groups.keys():
        # Loop over states in typ1
        for ind1, st1 in enumerate(groups[typ1]):

            # Loop over types
            for typ2 in types_tmp:
                # Loop over states in typ2
                for ind2, st2 in enumerate(groups[typ2]):

                    # If a ground state ...
                    if typ1 == 'grd':
                        # ... is interacting with a spectator state
                        if typ2 == 'ion' and ind2 == settings['spectator']['ion_ind']:
                            st2 = groups[typ2][settings['spectator']['ion_ind']]
                            for nn2 in range(
                                settings['spectator']['n_spec_min'], 
                                settings['spectator']['n_spec_max']+1
                                ):
                                _calc_ce(
                                    fac=fac, name = name,
                                    st1 = [st1],
                                    st2 = ([st2], nn2)
                                    )
                            # Highly excited state
                            if settings['spectator']['n_ext'] is not None:
                                _calc_ce(
                                    fac=fac, name=name,
                                    st1 = [st1],
                                    st2 = ([st2], settings['spectator']['n_ext'])
                                    )

                        # ... is interaction with a ground or excited state
                        elif typ2 in ['grd','exc']:
                            _calc_ce(fac=fac, name=name, st1=[st1], st2=[st2])

                    # If an excitation state ...
                    elif typ1 == 'exc':
                        # If skipping collisional configuration interactions
                        if 'skip_CI' in settings['restrictions'].keys():
                            continue

                        # ... is interacting with an excited state
                        if typ2 == 'exc':
                            # Considers just higher or equal quantum number
                            if ind1 <= ind2:
                                # Only calculate coll. exc. from at most the lowest excited states
                                if ind1 <= lim1 and ind2 <= lim2:
                                    _calc_ce(fac=fac, name=name, st1=[st1], st2=[st2])

        # Makes sure we don't improperly order the energy levels
        types_tmp.remove(typ1)

    # Print the time this took
    stop = time.time()
    s='Done %10.3E s'%(stop-start)
    print(s)

# Calculates transition rates
def _calc_ce(
    fac = None,
    name = None,
    st1 = None,
    st2 = None,
    ):

    # Calculates collisional excitation cross-sections
    print(str(st1)+' --> '+str(st2))
    fac.CETable(name+'b.ce',st1, st2)
    # Clears memory
    fac.Reinit(radial = 1, excitation = 1)