'''

_pt1.py is a module meant to utilize the first step of running FAC
to calculate the rates for various atomic processes as faciliated by the
"fac" module

NOTES:
    1) "pfac" -> python-version of FAC code base
    2) Based off of demo in FeL/fe.py, NeL/ne.py, https://www-amdis.iaea.org/FAC/
    3) This script is modeled off the provided pfac.atom.atomic_data() 
        specialized example
    4) This script is intended to calculate the atomic process of H-, He-, Li-like Ar
        to compare against FAC files shared to me by fsciortino

OUTPUTS:
    1) The "fac" module oututs files in the form SS##a.zz, SS##b.zz where
        i) SS is the atomic symbol
        ii) ## is the number of electrons
        iii) a -> ASCII format, b -> binary format
        iv) zz -> file type
            1) .en -> energy levels by fac.Structure
            2) .tr -> radiative transition rates by fac.TRTable
            3) .ce -> collisional excitation cross-sections by fac.CETable
            4) .ci -> collisional ionization cross-sections by fac.CITable
            5) .ai -> autoionization rates by fac.AITable
            6) .rr -. radiative recombination and photoionization cross-sections by fac.RRTable


TO DO:
    1) Generalize configuration method 
        i) Currently seems fine up to Ne-like
        ii) No innershell interactions right now, 
            just promoting highest ground quantum number
    2) Validate results against shared FAC data (Ar, Ca, W)

cjperks
03/13/23


'''

# Modules
import numpy as np
from pfac import fac
import time
import os

__all__ = [
    'rates'
    ]

##########################################################
#
#                     Main
#
##########################################################

def rates(
    # Input
    ion=None,           # Name of ion species
    save_path=None,     # path/to/save
    # Configuration Settings
    nele=1,             # Number of electrons
    settings = None,    # Dictionary of modeling settings
    #nmax=2,             # Maximum excitation quantum number
    #n_spec_min=10,      # Minimum spectator electron quantum number
    #n_spec_max=11,      # Maximum spectator electron quantum number
    # Physics Contols
    physics=None,       # If None -> include all
    ):

    # Sets the standard values for mass, charge, nuclear radius
    fac.SetAtom(ion)

    # Save file name
    name = os.path.join(
        save_path,
        '%s%02d'%(ion, nele)
        )

    # Physics to include
    if physics is None:
        physics = ['en', 'tr', 'ce', 'rr', 'ci', 'ai']

    # Calculates energy levels
    groups = _states(
        name=name,
        nele=nele,
        #nmax=nmax,
        settings=settings,
        )

    # Calculates energy levels
    if 'en' in physics:
        _en(
            name=name,
            groups=groups,
            n_spec_min=n_spec_min,
            n_spec_max=n_spec_max,
            )

    # Calculates radiative transition rates
    if 'tr' in physics:
        _tr(
            name=name,
            groups=groups,
            n_spec_min=n_spec_min,
            n_spec_max=n_spec_max,
            )

    # Calculates collisional excitation cross-sections
    if 'ce' in physics:
        _ce(
            name=name,
            groups=groups,
            n_spec_min=n_spec_min,
            n_spec_max=n_spec_max,
            )

    # Calculates radiative recombination rates
    if 'rr' in physics:
        _rr(
            name=name,
            groups=groups,
            n_spec_min=n_spec_min,
            n_spec_max=n_spec_max,
            )

    # Calculates collisional ionization cross-section
    if 'ci' in physics:
        _ci(
            name=name,
            groups=groups,
            n_spec_min=n_spec_min,
            n_spec_max=n_spec_max,
            )

    # Calculates  autoionization rates
    if 'ai' in physics:
        _ai(
            name=name,
            groups=groups,
            n_spec_min=n_spec_min,
            n_spec_max=n_spec_max,
            )

    # Builds an energy level table in memory from file
    fac.MemENTable(name+'b.en')

    # Calculates Maxwellian-averaged rate coefficients
    #temps = list(np.logspace(np.log10(8.62e0), np.log10(8.62e4), 41))
    #fac.MaxwellRate(name+'b.ce', name+'ce.mr', -1, -1, temp)
    #fac.MaxwellRate(name+'b.rr', name+'rr.mr', -1, -1, temp)
    #fac.MaxwellRate(name+'b.ci', name+'ci.mr', -1, -1, temp)

    # Converts binary file to ASCII
    #   Optional arg: v=1 -> verbose
    fac.PrintTable(name+'b.en', name+'a.en', 1)
    for ph in physics:
        fac.PrintTable(name+'b.'+ph, name+'a.'+ph, 1)

    # Clears all memory
    fac.Reinit(0)
        
##########################################################
#
#                     Utilities
#
##########################################################

# States of interest
def _states(
    name=None,
    nele=None,
    #nmax=None,
    settings=None,
    ):

    print('Generating states of interest')

    groups = {
        'grd': [],
        'exc': [],
        'ion': [],
        }

    # Loop over group types
    for typ in groups.keys():
        print(typ)
        # If n=1 shell is unfilled
        if nele < 3:
            # Filling the n=1 shell
            struct = '1*%d'%(nele)
            nmin = 2

        # If n=2 shell is unfilled
        elif nele >=3 and nele <=10:
            # Filling the n=2 shell
            struct = '1*2 2*%d'%(nele-2)
            nmin = 3

        # If considering a ground state
        if typ == 'grd':
            # Defines ground state
            groups[typ].append('grd.1')
            #fac.Config(groups[typ][-1], struct)
            print(struct)

        # If considering an excited state
        elif typ == 'exc':
            # --- Singly excited states --- #
            # Loop over excitation quantum numbers
            for nn in range(nmin,settings['single']['n_max']+1):
                if nele == 1 or nele == 3:
                    # Bringing an electron from ground to the n=nn shell
                    tmp = struct[:-3]

                else:
                    # Bringing an electron from ground to the n=nn shell
                    tmp = struct[:-1]
                    tmp += str(int(struct[-1])-1)
                    tmp += ' '

                # Excites an electron
                tmp += '%d*1'%(nn)

                # Defines this excited state
                groups[typ].append('exc.%d'%(nn))
                #fac.Config(groups[typ][-1], tmp)
                print(tmp)

            # --- Doubly excited states --- #
            if 'double' in settings.keys() and nele >1:
                # Loop over the lower excited electron
                for nn1 in range(
                    settings['double']['n1_min'],
                    settings['double']['n1_max']+1
                    ):
                    # Loop over the higher excited electron
                    for nn2 in range(
                        settings['double']['n2_min'],
                        settings['double']['n2_max']+1
                        ):

                        # If looking at innershell excitation
                        if nn1 <= nmin:
                            bb = [int(xx[0]) for xx in struct.split(' ')]
                            ee = [int(xx[-1]) for xx in struct.split(' ')]
                            ind = np.where(nn1 in bb)[0][0]

                            # Promotes innershell electron
                            tmp = ''
                            for ii,nn in enumerate(range(1,nmin)):
                                if ii < ind:
                                    tmp += str(bb[ii])+'*'+str(ee[ii])+' '

                                elif ii == ind:
                                    tmp += str(bb[ii])+'*'+str(ee[ii]-1)+' '

                                elif ii > ind:
                                    tmp += str(bb[ii])+'*'+str(ee[ii]+1)+' '
                            tmp = tmp[:-1]

                        else:
                            if nele == 1 or nele == 3:
                                # Bringing an electron from ground to the n=nn shell
                                tmp = struct[:-4]

                            else:
                                # Bringing an electron from ground to the n=nn shell
                                tmp = struct[:-1]
                                tmp += str(int(struct[-1])-1)
                                tmp += ' '
                            tmp += str(nn1)+'*1'

                        # Adds second excited electron
                        if nn2 == nn1:
                            tmp2 = tmp
                        elif nn2 > nn1:
                            tmp2 = tmp[:-1]
                            tmp2 += str(int(tmp[-1])-1)
                            tmp2 += ' %i*1'%(nn2)

                        # Error check
                        elif nn2 < nn1:
                            continue

                        # Defines this excited state
                        groups[typ].append('exc.%d.%d'%(nn1, nn2))
                        #fac.Config(groups[typ][-1], tmp2)
                        print(tmp2)

        # If considering an ionized state
        elif typ == 'ion':
            if nele <= 3:
                ni_min = 1
            elif nele <=11:
                ni_min = 2

            if nele == 1 or nele == 3:
                # Removing an electron from ground
                tmp = struct[:-4]

            else:
                # Removing an electron from ground
                tmp = struct[:-1]
                tmp += str(int(struct[-1])-1)

            # --- Outershell excitation --- #
            for nni in range(ni_min, settings['ioniz']['n_outer_max']+1):
                if nni == int(ni_min):
                    tmpi = tmp

                else:
                    if int(tmp[-1])-1 >0:
                        tmpi = tmp[:-1]
                        tmpi += str(int(tmp[-1])-1)
                    else:
                        tmpi = tmp[:-4]
                    tmpi += ' %i*1'%(nni)

                # Defines this ionized state
                print('outer')
                groups[typ].append('ion.%d'%(nni))
                #fac.Config(groups[typ][-1], tmp)
                print(tmpi)

            # --- Innershell excitation --- #
            for inn in settings['ioniz']['n_inner']:
                bb = [int(xx[0]) for xx in tmp.split(' ')]
                ee = [int(xx[-1]) for xx in tmp.split(' ')]
                ind = np.where(int(inn) in bb)[0][0]

                for nni in range(inn+1, settings['ioniz']['n_inner_max']+1):
                    print(nni)
                    tmpi = ''
                    for ii in range(len(ee)):
                        if ii < ind:
                            tmpi += str(bb[ii])+'*'+str(ee[ii])+' '
                        elif ii == ind:
                                tmpi += str(bb[ii])+'*'+str(ee[ii]-1)+' '

                        elif ii > ind:
                            if nni <= ee[-1]:
                                tmpi += str(bb[ii])+'*'+str(ee[ii]+1)+' '
                            else:
                                tmpi += str(bb[ii])+'*'+str(ee[ii])+' '


                    # Defines this ionized state
                    groups[typ].append('ion.%d'%(nni))
                    #fac.Config(groups[typ][-1], tmp)
                    print('inner')
                    print(tmpi)
                '''
                if nele == 1 or nele == 3:
                    # Removing an electron from ground
                    tmp = struct[:-4]

                else:
                    # Removing an electron from ground
                    tmp = struct[:-1]
                    tmp += str(int(struct[-1])-1)

                for nni in range(int(ni_min), settings['ioniz']['n_max']+1):
                    if nni == int(ni_min):
                        tmpi = tmp

                    else:
                    if nele == 2 or nele == 4:
                        bb = [int(xx[0]) for xx in tmp.split(' ')]
                        ee = [int(xx[-1]) for xx in tmp.split(' ')]
                        ind = np.where(int(nni-1) in bb)[0][0]

                        # Promotes innershell electron
                        tmpi = ''
                        for ii,nn in enumerate(range(1,int(nni))):
                            if ii < ind:
                                tmpi += str(bb[ii])+'*'+str(ee[ii])+' '

                            elif ii == ind:
                                tmpi += str(bb[ii])+'*'+str(ee[ii]-1)+' '

                            elif ii > ind:
                                tmpi += str(bb[ii])+'*'+str(ee[ii]+1)+' '
                        tmpi = tmpi[:-1]
                    else:
                        tmpi = tmp[:-1]
                        tmpi += str(int(tmp[-1])-1)
                        tmpi += ' %i*1'%(int(nni+1))
                '''


    # Stores configurations
    #fac.ListConfig()
    #fac.ListConfig(name+'a.cfg')
    
    # Optimizes the radial potential
    #   Purpose of this routine is to remove some error in energy level calcs
    #   due to using a single central potential for all configurations
    #fac.ConfigEnergy(0)         # Calculates the average energy of configurations
    #fac.OptimizeRadial(['grd.1']) # Obtaines optimal radial potential for list of cofigs
    #fac.ConfigEnergy(1)         # Recalculates average energy after optimized potential

    # Output
    return groups
                  
# Energy levels
def _en(
    name=None,
    groups=None,
    n_spec_min=None,
    n_spec_max=None,
    ):

    # Track time spent
    start = time.time()
    start0 = start

    print('EN: Calc energy levels')

    # Loop over groups
    for typ in groups.keys():
        # If considering a ground states
        if typ == 'grd':
            print('Structure: ground complex')
            print([groups[typ][-1], groups['exc'][0]])
            # Calculates atomic structure with first excited state
            fac.Structure(name+'b.en', 
                [groups[typ][-1], groups['exc'][0]]
                )

            # Precalculates angular coefficients with first excited state
            fac.PrepAngular([groups[typ][-1], groups['exc'][0]])

        # If considering an ionized state
        elif typ == 'ion':
            print('Structure: ionized complex')
            print([groups[typ][-1]])
            # Calculates atmoic structure
            fac.Structure(name+'b.en', 
                [groups[typ][-1]]
                )

            # Precalculates angular coefficients
            fac.PrepAngular([groups[typ][-1]])
            fac.PrepAngular([groups[typ][-1]], 
                [groups['grd'][-1]]
                )

            # Spectator electron settings
            fac.SetRecPWOptions(n_spec_min)
            fac.SetRecSpectator(n_spec_min,n_spec_min-1)

        # If considering an excited state
        elif typ == 'exc':
            print('Structure: excited complexes')
            # Loop over excited states excluding the first
            for st in groups[typ][1:]:
                print([st])
                # Calculates atomic structure
                fac.Structure(name+'b.en', 
                    [st]
                    )

                # Precalculates angular coefficients
                fac.PrepAngular([groups['ion'][-1]], 
                    [st]
                    )
                fac.PrepAngular([groups['grd'][-1]], 
                    [st]
                    )


    # Constructs spectator electron for recombining states
    for st in groups['ion']:
        for nn in range(n_spec_min, n_spec_max+1):
            print(([st], nn))
            fac.RecStates(name+'b.en',
                [st],
                nn
                )

    # Print the time this took
    stop = time.time()
    s='Done %10.3E s'%(stop-start)
    print(s)

# Radiative transition rates
def _tr(
    name=None,
    groups=None,
    n_spec_min=None,
    n_spec_max=None,
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
                        # ... is interacting with an ionized state
                        if typ2 == 'ion':
                            # Loop over spectator electrons
                            for nn2 in range(n_spec_min, n_spec_max+1):
                                # Calculates radiative transition rates
                                print(str([st1])+' --> '+str(([st2],nn2)))
                                fac.TRTable(name+'b.tr',
                                    [st1],
                                    ([st2], nn2)
                                    )
                                # Clears memory
                                fac.Reinit(radial = 1, excitation = 1)

                        # ... is interaction with a ground or excited state
                        else:
                            # Calculates radiative transition rates
                            print(str([st1])+' --> '+str([st2]))
                            fac.TRTable(name+'b.tr',
                                [st1],
                                [st2]
                                )
                            # Clears memory
                            fac.Reinit(radial = 1, excitation = 1)

                    # If an excitation state ...
                    elif typ1 == 'exc':
                        # ... is interacting with an excited state
                        if typ2 == 'exc':
                            # Considers just higher or equal quantum number
                            if ind1 <= ind2:
                                # Calculates radiative transition rates
                                print(str([st1])+' --> '+str([st2]))
                                fac.TRTable(name+'b.tr',
                                    [st1],
                                    [st2]
                                    )
                                # Clears memory
                                fac.Reinit(radial = 1, excitation = 1)

                        # ... is interacting with an ionized state
                        else:
                            # Just consider first excited state interacting with ionized state
                            if ind1 == 0:
                                # Loop over spectator electrons
                                for nn2 in range(n_spec_min, n_spec_max+1):
                                    # Calculates radiative transition rates
                                    print(str([st1])+' --> '+str(([st2],nn2)))
                                    fac.TRTable(name+'b.tr',
                                        [st1],
                                        ([st2], nn2),
                                        -1)
                                    # Clears memory
                                    fac.Reinit(radial = 1, excitation = 1)
                            
                    # If an ionized state ...
                    else:
                        # ... is interacting with an ionized state
                        for nn1 in range(n_spec_min, n_spec_max+1):
                            for nn2 in range(nn1, n_spec_max+1):
                                # Calculates radiative transition rates
                                print(str(([st1],nn1))+' --> '+str(([st2],nn2)))
                                fac.TRTable(name+'b.tr',
                                    ([st1], nn1),
                                    ([st2], nn2),
                                    -1)
                                # Clears memory
                                fac.Reinit(radial = 1, excitation = 1)

        # Makes sure we don't improperly order the energy levels
        types_tmp.remove(typ1)

    # Print the time this took
    stop = time.time()
    s='Done %10.3E s'%(stop-start)
    print(s)

# Collisional excitation cross-sections
def _ce(
    name=None,
    groups=None,
    n_spec_min=None,
    n_spec_max=None,
    ):

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
                        # ... is interacting with an ionized state
                        if typ2 == 'ion':
                            # Loop over spectator electrons
                            for nn2 in range(n_spec_min, n_spec_max+1):
                                # Consider only the lowest energy spectator electron
                                if nn2 == n_spec_min:
                                    # Calculates collisional excitation cross-sections
                                    print(str([st1])+' --> '+str(([st2],nn2)))
                                    fac.CETable(name+'b.ce',
                                        [st1],
                                        ([st2], nn2)
                                        )
                                    # Clears memory
                                    fac.Reinit(radial = 1, excitation = 1)

                        # ... is interaction with a ground or excited state
                        else:
                            # Calculates collisional excitation cross-sections
                            print(str([st1])+' --> '+str([st2]))
                            fac.CETable(name+'b.ce',
                                [st1],
                                [st2]
                                )
                            # Clears memory
                            fac.Reinit(radial = 1, excitation = 1)

                    # If an excitation state ...
                    elif typ1 == 'exc':
                        # ... is interacting with an excited state
                        if typ2 == 'exc':
                            # Considers just higher or equal quantum number
                            if ind1 <= ind2:
                                # Only calculate coll. exc. from at most the first excited state
                                if ind1 == 0:
                                    # Calculates collisional excitation cross-sections
                                    print(str([st1])+' --> '+str([st2]))
                                    fac.CETable(name+'b.ce',
                                        [st1],
                                        [st2]
                                        )
                                    # Clears memory
                                    fac.Reinit(radial = 1, excitation = 1)

        # Makes sure we don't improperly order the energy levels
        types_tmp.remove(typ1)

    # Print the time this took
    stop = time.time()
    s='Done %10.3E s'%(stop-start)
    print(s)

# Radiative recombination rates
def _rr(
    name=None,
    groups=None,
    n_spec_min=None,
    n_spec_max=None,
    ):

    # Track time spent
    start = time.time()
    start0 = start

    print('RR: Calc radiative recombination')

    # Loop over types
    for typ1 in groups.keys():
        # Loop over states
        for st1 in groups[typ1]:

            # Loop over ionized states
            for st2 in groups['ion']:

                if typ1 == 'ion':
                    # Loop over spectator electrons
                    for nn1 in range(n_spec_min, n_spec_max+1):
                        # Calculates rad. recomb and photoionization rates
                        print(str(([st1],nn1))+' --> '+str([st2]))
                        fac.RRTable(name+'b.rr',
                            ([st1], nn1),
                            [st2]
                            )
                        fac.SetPEGrid(0)
                        fac.SetUsrPEGrid(0)
                        fac.SetRRTEGrid(0)

                        # Clears memory
                        fac.Reinit(radial = 1, recombination = 1, ionization = 1)

                else:
                    # Calculates rad. recomb and photoionization rates
                    print(str([st1])+' --> '+str([st2]))
                    fac.RRTable(name+'b.rr',
                        [st1],
                        [st2]
                        )
                    fac.SetPEGrid(0)
                    fac.SetUsrPEGrid(0)
                    fac.SetRRTEGrid(0)

                    # Clears memory
                    fac.Reinit(radial = 1, recombination = 1, ionization = 1)

    # Print the time this took
    stop = time.time()
    s='Done %10.3E s'%(stop-start)
    print(s)

# Collisional ionization cross-section
def _ci(
    name=None,
    groups=None,
    n_spec_min=None,
    n_spec_max=None,
    ):

    # Track time spent
    start = time.time()
    start0 = start

    print('CI: Calc collisional ionization')

    # Loop over types
    for typ1 in groups.keys():
        # Loop over states
        for st1 in groups[typ1]:

            # Loop over ionized states
            for st2 in groups['ion']:

                if typ1 == 'ion':
                    # Loop over spectator electrons
                    for nn1 in range(n_spec_min, n_spec_max+1):
                        # Calculates collisional excitation cross-sections
                        print(str(([st1],nn1))+' --> '+str([st2]))
                        fac.CITable(name+'b.ci',
                            ([st1], nn1),
                            [st2]
                            )

                        # Clears memory
                        fac.Reinit(radial = 1, recombination = 1, ionization = 1)

                else:
                    # Calculates collisional excitation cross-sections
                    print(str([st1])+' --> '+str([st2]))
                    fac.CITable(name+'b.ci',
                        [st1],
                        [st2]
                        )

                    # Clears memory
                    fac.Reinit(radial = 1, recombination = 1, ionization = 1)

    # Print the time this took
    stop = time.time()
    s='Done %10.3E s'%(stop-start)
    print(s)

# Autoionization rates
def _ai(
    name=None,
    groups=None,
    n_spec_min=None,
    n_spec_max=None,
    ):

    # Track time spent
    start = time.time()
    start0 = start

    print('AI: Calc autoionization')

    # Loop over types
    for typ1 in list(groups.keys()):
        if typ1 == 'grd':
            continue
        # Loop over states
        for st1 in groups[typ1]:

            # Loop over ionized states
            for st2 in groups['ion']:

                if typ1 == 'ion':
                    # Loop over spectator electrons
                    for nn1 in range(n_spec_min, n_spec_max+1):
                        # Calculates autoionization rates
                        print(str(([st1],nn1))+' --> '+str([st2]))
                        fac.AITable(name+'b.ai',
                            ([st1], nn1),
                            [st2]
                            )

                        # Clears memory
                        fac.Reinit(radial = 1, recombination = 1)

                else:
                    # Calculates autoionization rates
                    print(str([st1])+' --> '+str([st2]))
                    fac.AITable(name+'b.ai',
                        [st1],
                        [st2]
                        )

                    # Clears memory
                    fac.Reinit(radial = 1, recombination = 1)

    # Print the time this took
    stop = time.time()
    s='Done %10.3E s'%(stop-start)
    print(s)


