'''

_ai_dr.py is a script meant to facilitate calculating
autoionization and dielectronic recombination using FAC

cjperks
May 9, 2024

'''

# Modules
import numpy as np
import sys, os
import copy
import time

__all__ = [
    '_run_ai'
    ]


######################################################
#
#        Autoionization/dielectronic recombination
#
######################################################

# Autoionization rates
def _run_ai(
    name = None,
    groups = None,
    nele = None,
    fac = None,
    ):
    '''
    NOTE: Assumes spectator states built as singly ionized
    '''

    # Ignore H-like
    if nele == 1:
        return

    # Track time spent
    start = time.time()
    start0 = start

    print('AI: Calc autoionization')

    # Loop over excited states
    for ind1, st1 in enumerate(groups['exc']):

        # Loop over ionized states
        for st2 in groups['ion']:

            # Ionizing from excitated state
            _calc_ai(fac=fac, name=name, st1=[st1], st2=[st2])

    # Print the time this took
    stop = time.time()
    s='Done %10.3E s'%(stop-start)
    print(s)

# Calculates autoionization rates
def _calc_ai(
    fac = None,
    name = None,
    st1 = None,
    st2 = None,
    ):

    # Calculates autoionization rates
    print(str(st1)+' --> '+str(st2))
    fac.AITable(name+'b.ai', st1, st2)

    # Clears memory
    fac.Reinit(radial = 1, recombination = 1)