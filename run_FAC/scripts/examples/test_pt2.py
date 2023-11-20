'''

test_pt2.py is a script meant to utilize the second step of running FAC
to run its collisional-radiative model ("crm" module) to calculate line 
intensities from atomic data calculated by "fac"

NOTES:
    1) "pfac" -> python-version of FAC code base
    2) Demo in FeL/spec.py, NeL/spec.py, https://www-amdis.iaea.org/FAC/

TO DO:
    1) Couple with Aurora impurity transport modeling -> nz_cm3



cjperks
03/13/23


'''

# Modules
import numpy as np
import matplotlib.pyplot as plt
from pfac import crm, fac

####################################################################
#
#           Defining plasma of interest
#
####################################################################

# Ion of interest
ion = 'Ar'

# Number of electrons of interest
nele = 1 # H-like
#nele = 2 # He-like
#nele = 3 # Li-like


# Electron kinetic profiles
Te_eV = 10e3 # [eV]
ne_cm3 = 1e14 # [cm^-3]

# Impurity kinetic profiles
ni_cm3 = 1e10 # [cm^-3], total impurity density
fz = 0.8 # charge state fractional abundance

# String of the ion charge state of interest
#   Used as file name
p = '%s%02d'%(ion, nele)


####################################################################
#
#           Collisional-radiative modeling
#
####################################################################

# crm.AddIon(nele, nz, 'file_prefix')
#   Add an ion to the spectral modeling
#       If want more than one, must add ion inc nele order
#   file_prefix -> base binary file name for 'fac' data files

# crm.SetBlocks([den, 'file_prefix'])
#   Reads the energy levels and sets superlevel blocks for spectral modeling
#   Optional:
#       den -> abundance of ion with one less eelctron than lowest charge state
#           added with crm.AddIon()
#       If den < 0 then only processes levels within the same ion

# crm.SetEleDist(i,...)
#   Sets the electron energy distribution
#   i -> type of distribution
#       0 -> Maxwellian
#       1 -> Gaussian
#       2 -> MaxPower
#       3 -> Power Law
#       -1 -> File
#   Other args are necessary params per the dist

# crm.SetIteration(a, [s,m])
#   Controls for population iterations
#   a -> accuracy (default 1e-4)
#   s -> stabilizer (default 0.75)
#   m -> max # of iterations (default 256)

# crm.SpecTable('file_name', [,rrc])
#   Outputs level populations and line emissivities
#   If rrc =1 -> incl radiative recomibnation continuum strength
#          =-1 -> only level populations

# crm.DumpRates('file_name', nele, type, [, imax [,a]])
#   Dump rate coefficients to a file
#   m -> which rate is outputted
#       =0 -> dump level indices, 2J values, energies, etc.
#       =1 -> dump radiative transition rates 
#       =2 -> dump two-photon transitions
#       =3 -> dump collisional excitation rates
#       =4 -> dump radiative recombination rates
#       =5 -> dump autoionization rates
#       =6 -> dump collisional ionization rates
#       transition rates in [1/s], rate coefficients in [1e-10 cm^3/s]
#   imax -> maximum levels to be included
#       <0 -> all
#   a =0 -> binary format, =1 -> ASCII format

# crm.PlotSpec('read_file', 'write_file, nele, t, e0, e1, de)
#   Prints spectrum of lines to 'write_file'
#   t -> line type
#       =0 -> all
#   e0 -> min energy of spectral range, [eV]
#   e1 -> max energy of spectral range, [eV]
#   de -> Gaussian FWHM, [eV]


# add the ion with nele number of electrons
crm.AddIon(nele, ni_cm3, '%sb'%p)

# Reads the data files
crm.SetBlocks(-1)

# Sets the electron distribution function
crm.SetEleDist(0, temp, -1, -1)

# Initializes radiative transition rates
#   If arg =1 then include photo-excitation (inverse process)
crm.SetTRRates(0)

# Initializes collisional excitation rate coefficients.
#   If arg =1 then include collisional de-excitation
crm.SetCERates(1)

# Initializes collisional ionization rate
#   If arg =1 then incl. 3-body recomb (not implemented yet)
crm.SetCIRates(0)

# abundance of the ion with nele number of electrons
crm.SetAbund(nele, fz)
crm.SetEleDensity(ne_cm3)

# Initializes superlevel blocks of the spectral model
crm.InitBlocks()

# Iteration controls
crm.SetIteration(1e-6, 0.5)

# Solver
crm.LevelPopulation()

# Prints output data file
crm.SpecTable(p+'b.sp', 0)

# Converts binary file to ASCII
#   Optional arg: v=1 -> verbose
crm.PrintTable(p+'b.sp', p+'a.sp')

# Saves collisional excitation rate coefficient data
crm.DumpRates(p+'a.rt3', nele, 3, -1, 1)

# Saves radiative transition rate coefficient data
crm.DumpRates(p+'a.rt1', nele, 1, -1, 1)

# Prints the Doppler broadened spectrum
#crm.PlotSpec(p+'b.sp', p+'spec', nele, 0, 3.7, 3.8, 4e3)