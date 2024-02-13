'''

This script is meant for learning how to use ChiantiPy
to model free-free, free-bound, and 2-photon continuum
emission

cjperks
jan 26, 204

'''

# Modules
import ChiantiPy.core as ch

# Useful constants
eV2K = 11604.5

# Plasma
Te = [3e3*eV2K] # [K]
ne = 1e14 # [1/cm3]

# emission measure
em = [1]

# wavelength grid
lamb = np.linspace(0.5, 10, 1000)

# Defines ion
## NOTE: spectroscopic notation
myIon = 'fe_25' # He-like Fe

# Init module
c = ch.continuum(
    myIon,
    temperature = Te,
    em = em,
    )

# Calculates free-free emission
c.freeFree(lamb)

# Calculates free-bound emission
c.freeBound(lamb)

# Calculates 2-photon
fe25=ch.ion(
    myIon, 
    Te, 
    ne, 
    em=em)
fe25.twoPhoton(lamb)

# Total emission
itemp = 0
total = (
    c.FreeFree['intensity'] 
    + c.FreeBound['intensity']
    + fe25.TwoPhoton['intensity']
    )

plt.figure()
plt.plot(lamb, c.FreeFree['intensity'],label='ff')
plt.plot(lamb, c.FreeBound['intensity'],label='fb')
plt.plot(lamb,fe25.TwoPhoton['intensity'],label='2 photon')
plt.plot(lamb, total, 'k', label='total')
plt.xlabel(c.FreeFree['xlabel'], fontsize=14)
plt.ylabel(c.FreeFree['ylabel'], fontsize=14)
plt.legend(loc='upper right', fontsize=14)
plt.title(' %s  T = %10.2e'%(fe25.IonStr, Te[itemp]), fontsize=14)
plt.ylim(bottom=0.)
plt.xlim([0., lamb[-1]])
plt.tight_layout