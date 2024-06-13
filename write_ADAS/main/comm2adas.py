'''

Script to manage reformatting a common data structure to writing
ADAS-format data

cjperks
June 13, 2024

'''

# Modules
import numpy as np
import os
from datetime import date

__all__ = [
    'write'
]

##########################################
#
#            Main
#
##########################################

# Formats data to ADAS-format
def write(
    # Atom
    sym = None,
    # File management
    filepath = None,
    filetype = None,
    source = None,
    # Data
    Te_eV = None, # [eV], dim(Te,)
    ne_cm3 = None, # [1/cm3], dim(ne,)
    data = None, # [cm3/s], dim(ne, Te, cs)
    ):

    # Get atomic information
    datom = _get_atom(sym=sym)

    # Init file write
    f = open(
        os.path.join(
            filepath,
            (
                filetype.lower() 
                + '_' +source.lower()
                + '_cjperks_' 
                + sym.lower() 
                + '.dat'
                )
            ),
        'w'
        )

    # Writes ACD/SCD
    if filetype in ['ACD', 'SCD']:
        _write_cs(
            f=f,
            filetype=filetype,
            Te_eV=Te_eV,
            ne_cm3=ne_cm3,
            data=data,
            datom=datom,
            source=source,
            )


# Writes effective ionization/recombination data file
def _write_cs(
    f=None,
    filetype=None,
    source=None,
    Te_eV=None,
    ne_cm3=None,
    data=None,
    datom=None,
    ):

    # Writes header
    f.write(
        (
            ('%i'%(datom['Zimp'])).rjust(5, ' ')
            +('%i'%(len(ne_cm3))).rjust(5, ' ')
            +('%i'%(len(Te_eV))).rjust(5, ' ')
            +('%i'%(1)).rjust(5, ' ')
            +('%i'%(datom['Zimp'])).rjust(5, ' ')
            +'     /'+datom['name']
            +'\n'
            )
        )
    f.write(' ---------------------------------------------------\n')

    # Writes density block
    _write_neTe_block(
        f=f,
        data=ne_cm3
        )

    # Writes temperature block
    _write_neTe_block(
        f=f,
        data=Te_eV
        )

    # Data type
    if filetype == 'ACD':
        label = 'RECOM'
    elif filetype == 'SCD':
        label = 'IONIZ'

    # Loop over charge states
    for cs in range(1,datom['Zimp']+1):
        _write_cs_block(
            f=f,
            cs=cs,
            data=data[:,:,cs-1],
            label=label
            )

    # Footer
    today = date.today().strftime("%b-%d-%Y")
    f.write('C-----------------------------------------------------------------------\n')
    f.write('C\n')
    f.write((
        'C\t'
        +label +' rate coefficient data\n'
        ))
    f.write((
        'C\t'
        +'SOURCED FROM:\t'
        + source+'\n'
        ))
    f.write('C\n')
    f.write('C\tPRODUCER:\t Conor Perks, cjperks@psfc.mit.edu\n')
    f.write('C\tDATE:    \t'+ today +'\n')
    f.write('C\n')
    f.write('C-----------------------------------------------------------------------')

    # Closes file
    f.close()

##########################################
#
#            Utilities
#
##########################################

# Writes ionizatio/recombination data block per charge state
def _write_cs_block(
    f=None,
    cs=None,
    data=None, # dim(ne,Te)
    label=None
    ):

    # Write header
    today = date.today().strftime("%m.%d.%y")
    f.write((
        ' -----------------'
        +'/ IPRT= 1\t'
        +'/ IGRD= 1\t'
        +'/'+label+'\t'
        +'/ Z1= %i\t'%(cs)
        +'/ DATE= '+today
        +'\n'
        ))

    # Number of columns
    ncols = 8

    # Number of rows per density point
    nrows = int(np.ceil(data.shape[0]/ncols))

    # Loop over temperature points
    for tt in range(data.shape[1]):
        # Loop over rows
        for rr in range(nrows):
            row = ''
            # Loop over columns
            for cc in range(ncols):
                try:
                    row += ('%0.5f'%(
                        np.log10(data[cc + rr*ncols, tt])
                        )).rjust(10, ' ')
                except:
                    row += ''

            # Writes row
            row += '\n'
            f.write(row)


# Writes block for ne/Te
def _write_neTe_block(
    f=None,
    data=None,
    ):

    # Number of columns
    ncols = 8

    # Number of rows
    nrows = int(np.ceil(len(data)/ncols))

    # Loop over rows
    for rr in range(nrows):
        row = ''
        # Loop over columns
        for cc in range(ncols):
            try:
                row += ('%0.5f'%(
                    np.log10(data[cc+rr*ncols])
                    )).rjust(10, ' ')
            except:
                row += ''

        # Write rows
        row += '\n'
        f.write(row)

# Gets atomic information
def _get_atom(
    sym = None,
    ):

    datom = {
        'Ar':{
            'name': 'ARGON',
            'Zimp': 18
            },
        'Kr':{
            'name': 'KRYPTON',
            'Zimp': 36
            }
        }

    # Output
    return datom[sym]