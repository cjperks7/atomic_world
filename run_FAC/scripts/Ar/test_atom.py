
from pfac import atom

Ar = atom.ATOM(1,asym='Ar', dir='test_atom/')

Ar.set_atom('Ar', dir='test_atom/')

Ar.set_configs()

Ar.exc_complex[0].cgroup[0].name

Ar.ion_complex.cgroup[0].name

atom.atomic_data(1,'Ar', dir='test_atom/')