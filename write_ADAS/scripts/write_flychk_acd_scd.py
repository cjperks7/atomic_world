'''

Script to write FLYCHK data into ADAS-format (adf11) acd,scd files

cjperks
Jun 13, 2024

'''

# Modules
from atomic_world.write_ADAS.main import flychk2comm as F2C
from atomic_world.write_ADAS.main import comm2adas as C2A
import os
from transport_world.run_AURORA.calc_imp import _calc_imp_dens as cid

# Ion
sp = 'Ar'

# file management
file_read = os.path.join(
    '/home/cjperks/work',
    '2008_SPARCxray/FLYCHK_NLTE'
    )
file_write = os.path.join(
    '/home/cjperks',
    'test'
    )

# Reformats FLYCK data to common data format
Te_eV, ne_cm3, acd, scd = F2C.convert(
    filepath = file_read,
    filetype = 'ion_rates',
    sp=sp
    )

# Write common data format to ADAS-format
C2A.write(
    # Atom
    sym = sp,
    # File management
    filepath = file_write,
    filetype = 'ACD',
    source = 'FLYCHK',
    # Data
    Te_eV = Te_eV, # [eV], dim(Te,)
    ne_cm3 = ne_cm3, # [1/cm3], dim(ne,)
    data = acd, # [cm3/s], dim(ne, Te, cs)
    )
C2A.write(
    # Atom
    sym = sp,
    # File management
    filepath = file_write,
    filetype = 'SCD',
    source = 'FLYCHK',
    # Data
    Te_eV = Te_eV, # [eV], dim(Te,)
    ne_cm3 = ne_cm3, # [1/cm3], dim(ne,)
    data = scd, # [cm3/s], dim(ne, Te, cs)
    )

# Plots charge state distribution
fz = cid.get_ion_bal(
    sp=sp,
    Te_eV = np.logspace(np.log10(1e2), np.log10(1e4), 100),
    ne_cm3 = np.ones(100),
    files = {
        'acd':file_write+'/acd_cjperks_ar.dat',
        'scd':file_write+'/scd_cjperks_ar.dat'
        },
    plot = True
    )