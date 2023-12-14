'''

test_pt1.py is a script meant to utilize the first step of running FAC
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
import matplotlib.pyplot as plt
from pfac import fac
import time


####################################################################
#
#           Defining ion of interest
#
####################################################################

# Ion of interest
ion = 'Ar'

# Number of electrons of interest
nele = 1 # H-like
#nele = 2 # He-like
#nele = 3 # Li-like

# Maximum excitation quantum number
nmax = 2

# Spectator electron quantum number settings for ionized state
n_spec_min = 10
n_spec_max = 11

# String of the ion charge state of interest
#   Used as file name
p = 'test_ai/%s%02d'%(ion, nele)

# Sets the standard values for mass, charge, nuclear radius
fac.SetAtom(ion)


####################################################################
#
#          Defining energy levels of interest
#
####################################################################

# Notation used in the following
# X.Y.Z -> 
#   X -> n, quantum number of state of interest
#   Y -> l, orbital angular momentum of state of interest
#   Z -> j, total angular momentum of state of interest

# fac.Config('name', 'group')
#   Note: 'group' = 'X*Y' goes through all permutations of Y electron at n =X quantum number
#       Otherwise use spectroscopic symbols
#       Can use brackets about the spectroscopic symbol to denote j
#           +/- denotes net 1/2 electron spin
#       Can stimulate multiple conditions
#           i.e. 3*10; 3s>0; 3p>5
#               Generates 10 electrons at n=3
#               at least in 3s state and at least 6 in 3p state

print('Generating states of interest')

groups = {
    'grd': [],
    'exc': [],
    'ion': [],
}

# Loop over group types
for typ in groups.keys():
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
        fac.Config(groups[typ][-1], struct)

    # If considering an excited state
    elif typ == 'exc':
        # Loop over excitation quantum numbers
        for nn in range(nmin,nmax+1):
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
            fac.Config(groups[typ][-1], tmp)

    # If considering an ionized state
    elif typ == 'ion':
        if nele == 1 or nele == 3:
            # Removing an electron from ground
            tmp = struct[:-4]

        else:
            # Removing an electron from ground
            tmp = struct[:-1]
            tmp += str(int(struct[-1])-1)

        # Defines this ionized state
        groups[typ].append('ion.0')
        fac.Config(groups[typ][-1], tmp)
                
        

# Defines the ground state
#fac.Config('grd.1', '1*1') # Ground state 1s1

# Other energy levels of interest ****************************
#
#
#
#fac.Config('exc.0.2.0','2*1') # first excited state 2s1
#fac.Config('exc2.p.+', '2[p+]1') # first excited state 2p(3/2) -> 2p+
#fac.Config('exc2.p.-', '2[p-]1') # first excited state 2p(1/2) -> 2p-
#fac.Config('exc.2', '2*1')

# ionization state
#fac.Config('ion.0', '')

# Optimizes the radial potential
#   Purpose of this routine is to remove some error in energy level calcs
#   due to using a single central potential for all configurations
fac.ConfigEnergy(0)         # Calculates the average energy of configurations
fac.OptimizeRadial(['grd.1']) # Obtaines optimal radial potential for list of cofigs
fac.ConfigEnergy(1)         # Recalculates average energy after optimized potential


####################################################################
#
#          Calculating energy levels of interest
#
####################################################################

# Track time spent
start = time.time()
start0 = start


# fac.Structure('file name', ['groups'])
#   Diagonalizes the Hamiltonians for configus in 'groups'
#       Default is to treat interaction between groups exactly
#   from what I gather, need to do the ground an next excited state together
#       very similar states should also be grouped together
#       else, can just do one group

# fac.PrepAngular(['p', 'q'])
#   Pre-calculates angular coefficients between states 'p' and'q'
#       If 'q' is not present then just the states in list 'p'
#   This routine should primarily be used when atomic states are construced with
#       RecStates, where the angular coefficients between the base states 
#       are used many times.


# Settings for the spectator electron in recombining states
#   fac.SetRecPWOptions(lmax)
#       lmax -> maximum orbital angular momentum, default 12
#   fac.SetRecSpectator(nmax, nfrozen)
#       Settings for minimum qunatum number of spectator electrons
#       nmax -> States with n>nmax are set by fac.RecStates() function
#       nfrozen -> States with n>nfrozen are treated with frozen core approx
#       Default for both is 8

# fac.RecStates('file_name', 'state', 'quantum #')
#   Constructs recombining states by adding a spectator electron onto basis state
#   Should only be used if free electron is captured by empty orbital

print('EN: Calc energy levels')

# Loop over groups
for typ in groups.keys():
    # If considering a ground states
    if typ == 'grd':
        # Calculates atomic structure with first excited state
        fac.Structure(p+'b.en', 
            [groups[typ][-1], groups['exc'][0]]
            )

        # Precalculates angular coefficients with first excited state
        fac.PrepAngular([groups[typ][-1], groups['exc'][0]])

    # If considering an ionized state
    elif typ == 'ion':
        # Calculates atmoic structure
        fac.Structure(p+'b.en', 
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
        # Loop over excited states excluding the first
        for st in groups[typ][1:]:
            # Calculates atomic structure
            fac.Structure(p+'b.en', 
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
        fac.RecStates(p+'b.en',
            [st],
            nn
            )


# Calculates atomic structure
#fac.Structure(p+'b.en', ['grd.1', 'exc.2'])
#fac.PrepAngular(['grd.1', 'exc.2'])

#fac.Structure(p+'b.en', ['ion.0'])
#fac.PrepAngular(['ion.0'])
#fac.PrepAngular(['ion.0'],['grd.0'])
#fac.SetRecPWOptions(10)
#fac.SetRecSpectator(10,9)

# Constructs spectator electron for recombining stats
#fac.RecStates(p+'b.en',['ion.0.0'],10)
#fac.RecStates(p+'b.en',['ion.0.0'],11)


# Print the time this took
stop = time.time()
s='Done %10.3E s'%(stop-start)
print(s)



####################################################################
#
#          Calculating radiative transition rates &
#               collisional exciation cross-sections
#
####################################################################

# fac.TRTable('file name', 'low', 'up', ['multipole'])
#   Calculates oscillator strength and radiative transition rates 
#       from 'up' to 'low'
#   'multipole' -> list of multipole types
#       Default for 'multipoles'=0, so sum over all
#       1 = dipole
#       2 = quadrapole
#       (-) = electric type
#       (+) = magnetic type

# fac.CETable('file name', 'low', 'up')
#   Calculates collision strength 
#       from 'low' to 'up'
#   Note only need to do once for excite to ionized + spectator

# fac.Reint(m**key)
#   Reinitializes some or all subsystems of FAC package
#       I guess this is to cut down on memory?????????????????????????????
#   mkeys
#       config >=0 -> clear existing configs and couplings
#       recouple >=0 -> clear all coupling interactions
#       structure >=0 -> clear energy level table and angular coefficients 
#       radial -> radial wavefunction and integrals
#           =0 -> clear radial orbital, potential, and all radial integral tables
#           =1 -> clear only continuum orbitals
#       excitation >=0 -> clear collision energy grid and radial integrals
#       recombination -> photoionization, recombination, and autoionization
#           =0 -> clear energy grid, radial integrals, and recombined complex table
#           =1 -> clear only energy grid and radial integrals
#       ionization >=0 -> clear energy grid and radial integrals
#       dbase -> database handling
#           =0 -> clear energy table in memory and reinit all database headers
#           =t >0 -> clear header of database type t

print('TR: Calc radiative transtion')

# List of state types
types_tmp = list(groups.keys())

# Loop over types
for typ1 in groups.keys():
    # Loop over states in typ1
    for ind1, st1 in enumerate(groups[typ1]):
        print('First')
        print(typ1)
        print(st1)

        # Loop over types
        for typ2 in types_tmp:
            # Loop over states in typ2
            for ind2, st2 in enumerate(groups[typ2]):
                print('Second')
                print(typ2)
                print(st2)

                # If a ground state ...
                if typ1 == 'grd':
                    # ... is interacting with an ionized state
                    if typ2 == 'ion':
                        # Loop over spectator electrons
                        for nn2 in range(n_spec_min, n_spec_max+1):
                            # Calculates radiative transition rates
                            fac.TRTable(p+'b.tr',
                                [st1],
                                ([st2], nn2)
                                )
                            # Clears memory
                            fac.Reinit(radial = 1, excitation = 1)

                    # ... is interaction with a ground or excited state
                    else:
                        # Calculates radiative transition rates
                        fac.TRTable(p+'b.tr',
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
                            fac.TRTable(p+'b.tr',
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
                                fac.TRTable(p+'b.tr',
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
                            fac.TRTable(p+'b.tr',
                                ([st1], nn1),
                                ([st2], nn2),
                                -1)
                            # Clears memory
                            fac.Reinit(radial = 1, excitation = 1)

    # Makes sure we don't improperly order the energy levels
    types_tmp.remove(typ1)


print('CE: Calc collisional excitation')

# List of state types
types_tmp = list(groups.keys())

# Loop over types
for typ1 in groups.keys():
    # Loop over states in typ1
    for ind1, st1 in enumerate(groups[typ1]):
        print('First')
        print(typ1)
        print(st1)

        # Loop over types
        for typ2 in types_tmp:
            # Loop over states in typ2
            for ind2, st2 in enumerate(groups[typ2]):
                print('Second')
                print(typ2)
                print(st2)

                # If a ground state ...
                if typ1 == 'grd':
                    # ... is interacting with an ionized state
                    if typ2 == 'ion':
                        # Loop over spectator electrons
                        for nn2 in range(n_spec_min, n_spec_max+1):
                            # Consider only the lowest energy spectator electron
                            if nn2 == n_spec_min:
                                # Calculates collisional excitation cross-sections
                                fac.CETable(p+'b.ce',
                                    [st1],
                                    ([st2], nn2)
                                    )
                                # Clears memory
                                fac.Reinit(radial = 1, excitation = 1)

                    # ... is interaction with a ground or excited state
                    else:
                        # Calculates collisional excitation cross-sections
                        fac.CETable(p+'b.ce',
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
                                fac.CETable(p+'b.ce',
                                    [st1],
                                    [st2]
                                    )
                                # Clears memory
                                fac.Reinit(radial = 1, excitation = 1)

    # Makes sure we don't improperly order the energy levels
    types_tmp.remove(typ1)



# Calculates radiative transition rates
#   & collisional excitation cross-sections
#   then clears memory
#fac.TRTable(p+'b.tr', ['grd.1'], ['grd.1'])
#fac.CETable(p+'b.ce', ['grd.1'], ['grd.1']) # ground to everything
#fac.Reinit(radial = 1, excitation = 1)

#fac.TRTable(p+'b.tr', ['grd.1'], ['exc.2'])
#fac.CETable(p+'b.ce', ['grd.1'], ['exc.2'])
#fac.Reinit(radial = 1, excitation = 1)

#fac.TRTable(p+'b.tr', ['grd.1'], (['ion.0'], 10))
#fac.CETable(p+'b.ce', ['grd.1'], (['ion.0'], 10))
#fac.Reinit(radial = 1, excitation = 1)

#fac.TRTable(p+'b.tr', ['grd.1'], (['ion.0'], 11))
#fac.Reinit(radial = 1, excitation = 1)

#fac.TRTable(p+'b.tr', ['exc.2'], ['exc.2'])
#fac.CETable(p+'b.ce', ['exc.2'], ['exc.2']) # can do first excited to higher
#fac.Reinit(radial = 1, excitation = 1)

#fac.TRTable(p+'b.tr', ['exc.2'], (['ion.0'],10), -1)
#fac.Reinit(radial = 1, excitation = 1)

#fac.TRTable(p+'b.tr', ['exc.2'], (['ion.0'],11), -1)
#fac.Reinit(radial = 1, excitation = 1)

#fac.TRTable(p+'b.tr', (['ion.0'], 10), (['ion.0'], 10), -1) # higher n to lowest n
#fac.Reinit(radial = 1, excitation = 1)

#fac.TRTable(p+'b.tr', (['ion.0'], 10), (['ion.0'], 11), -1)
#fac.Reinit(radial = 1, excitation = 1)

#fac.TRTable(p+'b.tr', (['ion.0'], 11), (['ion.0'], 11), -1) # same n to same n
#fac.Reinit(radial = 1, excitation = 1)

fac.MaxwellRate(p+'b.ce', p+'ce.mr', -1, -1, [1370, 3e3, 10e3])

fac.MaxwellRate(p+'b.rr', p+'rr.mr', -1,-1, [1370, 3e3, 10e3])

fac.MaxwellRate(p+'b.ci', p+'ci.mr', -1,-1, [1370, 3e3, 10e3])

####################################################################
#
#          Calculating radiative recombination 
#               and photoionization cross-sections
#               & collisional ionization cross-sections
#
####################################################################

# fac.RRTable('file name', 'bound', 'free', 'multipole')
#   Calculates bound-free differential oscillator strengths

# fac.CITable('file name', 'bound', 'free')
#   Calculates the collision strength from ionizing 
#       'bound' group to the 'free'

# fac.SetPEGrid(#)
#   Sets the free electron energy grid for RR, PI, and AI calcs
#   System default settings = 0

# fac.SetUsrPEGrid(#)
#   Sets outputted free electron energy grid for RR and PI calcs
#   System default settings = 0

# fac.SetRRTEGrid(#)
#   Sets transition energy grid for RR and PI calcs
#       Not involved for electric dipole transitions
#   System default settings = 0

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
                    fac.RRTable(p+'b.rr',
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
                fac.RRTable(p+'b.rr',
                    [st1],
                    [st2]
                    )
                fac.SetPEGrid(0)
                fac.SetUsrPEGrid(0)
                fac.SetRRTEGrid(0)

                # Clears memory
                fac.Reinit(radial = 1, recombination = 1, ionization = 1)

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
                    fac.CITable(p+'b.ci',
                        ([st1], nn1),
                        [st2]
                        )

                    # Clears memory
                    fac.Reinit(radial = 1, recombination = 1, ionization = 1)

            else:
                # Calculates collisional excitation cross-sections
                fac.CITable(p+'b.ci',
                    [st1],
                    [st2]
                    )

                # Clears memory
                fac.Reinit(radial = 1, recombination = 1, ionization = 1)

# Calculate rad. recomb and photoion rates
#   & collisional excitation cross-sections
#   then clear memory
#fac.RRTable(p+'b.rr',['grd.1'],['ion.0'])
#fac.SetPEGrid(0)
#fac.SetUsrPEGrid(0)
#fac.SetRRTEGrid(0)
#fac.CITable(p+'b.ci', ['grd.1'], ['ion.0'])
#fac.Reinit(radial = 1, recombination = 1, ionization = 1)


#fac.RRTable(p+'b.rr',['exc.2'],['ion.0'])
#fac.SetPEGrid(0)
#fac.SetUsrPEGrid(0)
#fac.SetRRTEGrid(0)
#fac.CITable(p+'b.ci', ['exc.2'], ['ion.0'])
#fac.Reinit(radial = 1, recombination = 1, ionization = 1)


#fac.RRTable(p+'b.rr',(['ion.0'], 10),['ion.0'])
#fac.SetPEGrid(0)
#fac.SetUsrPEGrid(0)
#fac.SetRRTEGrid(0)
#fac.CITable(p+'b.ci', (['ion.0'], 10), ['ion.0'])
#fac.Reinit(radial = 1, recombination = 1, ionization = 1)


#fac.RRTable(p+'b.rr',(['ion.0'], 11),['ion.0'])
#fac.SetPEGrid(0)
#fac.SetUsrPEGrid(0)
#fac.SetRRTEGrid(0)
#fac.CITable(p+'b.ci', (['ion.0'], 11), ['ion.0'])
#fac.Reinit(radial = 1, recombination = 1, ionization = 1)



####################################################################
#
#          Calculating autoionization rates
#
####################################################################

# fac.AITable('file name', 'bound', 'free')
#   Calculates the autoionization rates between bound and free config

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
                    fac.AITable(p+'b.ai',
                        ([st1], nn1),
                        [st2]
                        )

                    # Clears memory
                    fac.Reinit(radial = 1, recombination = 1)

            else:
                # Calculates autoionization rates
                fac.AITable(p+'b.ai',
                    [st1],
                    [st2]
                    )

                # Clears memory
                fac.Reinit(radial = 1, recombination = 1)

#fac.AITable(p+'b.ai', ['exc.2'], ['ion.0'])
#fac.Reinit(radial = 1, recombination = 1)

#fac.AITable(p+'b.ai', (['ion.0'], 10), ['ion.0'])
#fac.Reinit(radial = 1, recombination = 1)

#fac.AITable(p+'b.ai', (['ion.0'], 11), ['ion.0'])
#fac.Reinit(radial = 1, recombination = 1)


####################################################################
#
#          Saving data in a readable format
#
####################################################################

# Builds an energy level table in memory from file
fac.MemENTable(p+'b.en')

# Converts binary file to ASCII
#   Optional arg: v=1 -> verbose
fac.PrintTable(p+'b.en', p+'a.en', 1)

fac.PrintTable(p+'b.tr', p+'a.tr', 1)

fac.PrintTable(p+'b.ce', p+'a.ce', 1)

fac.PrintTable(p+'b.rr', p+'a.rr', 1)

fac.PrintTable(p+'b.ci', p+'a.ci', 1)

fac.PrintTable(p+'b.ai', p+'a.ai', 1)

# Clears all memory
fac.Reinit(0)