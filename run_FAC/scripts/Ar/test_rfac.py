

from pfac import rfac

en = rfac.read_en('test/Ar01a.en')

# Dictionaryies of level types
en[1] # dim(types,) -> (grd/exc, ion, spectators)

# level index
en[1][xx].keys()
en[1][xx]['ILEV']

# level energies
en[1][xx]['ENERGY']

tr = rfac.read_tr('test/Ar01a.tr')

# Dictionary of interaction types
tr[1] # dim(types,) -> 

# Einstein coefficients
tr[1][0]['rate'] # [1/s]

ce = rfac.read_ce('test/Ar01a.ce')

rr = rfac.read_rr('test/Ar01a.rr')

ci = rfac.read_ci('test/Ar01a.ci')

# python object version
dd = rfac.FLEV('test/Ar01a.en')

dd.e # -> energy
dd.v # -> VNL
dd.s # -> LS name
dd.n # -> JJ name
dd.e0 # -> ground energy
dd.ei # -> ionization energy
dd.ib # -> base state
dd.nele # -> number electrons


# --- Spectral modeling data --- #

# Spectral line data
sp = rfac.read_sp('test/Ar01a.sp')


# Line rate data
rt1 = rfac.read_rt('test/Ar01a.rt1')

