'''

_write_PEC.py is a module meant to manage writing adf15-format
PEC data files using FAC atomic data

cjperks
Dec 14, 2023

TO DO:
    1) How to get 2S+1 values?
        -> How also do you differentiate in FAC?

'''

# Modules
from atomic_world.run_FAC.main import _crm as crm
import os
import scipy.constants as cnt

__all__ = [
    'write_PEC'
    ]

######################################################
#
#                   Main
#
######################################################

def write_PEC(
    # Ion species
    ion = None,
    ZZ = None,
    nz_cm3 = None,
    # Plasma conditions
    Te_array = None, # [eV]
    ne_array = None, # [1/cm3]
    # Simulation settings
    path_fac = None,
    nsim = None,
    # File writing controls
    lamb_min = None, # [A]
    lamb_max = None, # [A]
    ):

    # Initialize dictionary to store data
    cpls = None     # dim(nd,nt)

    # Energy -> wavelength
    hc = cnt.h*cnt.c/cnt.e*1e10 # [eV*AA]

    # Loop over temperature points
    for nt in np.arange(len(Te_array)):
        # Formats FAC data
        data, lvls = crm.prep(
            ion = ion,
            # FAC controls
            path_fac = path_fac,
            nele = len(nz_cm3) -1 -ZZ,
            nsim = nsim,  # Maximum quantum number to solve for population
            mode = 'mr',
            Te_eV = Te_array[nt],
            )

        # Loop over density points
        for nd in np.arange(len(ne_array)):
            # Defines wants transitions to get PECs for
            if cpls is None:
                cpls = {}

                # Loop over upper levels to consider
                for upr in data.keys():
                    if data[upr]['E_eV'] > 0:
                        cpls[upr] = {}
                    else:
                        continue

                    # Loop over lower levels to consider
                    for lwr in lvls:
                        # Determines if transition is within wavelength range of interest
                        lamb = hc/(data[upr]['E_eV'] - data[lwr]['E_eV'])
                        if ((lamb >= lamb_min)
                            and (lamb <= lamb_max)
                            ):
                            cpls[upr][lwr] = {}
                            cpls[upr][lwr]['lambda_A'] = lamb

                            # PEC types to include
                            cpls[upr][lwr]['EXCIT'] = np.zeros((len(ne_array), len(Te_array)))
                            cpls[upr][lwr]['RECOM'] = np.zeros((len(ne_array), len(Te_array)))
                            cpls[upr][lwr]['IONIZ'] = np.zeros((len(ne_array), len(Te_array)))

            # Solves relative population balance
            pops, mat, rhs = crm.pops(
                data = data,
                lvls = lvls,    # [list], levels to solve for
                nsim = nsim,    
                # Plasma
                ne_cm3 = ne_array[nd],
                nz_cm3 = fz[0,:],  # dim(cs,)
                ZZ = ZZ,
                ) # [1/cm3]

            # Stores PECs
            ntrans = 0
            for upr in cpls.keys():
                for lwr in cpls[upr].keys():
                    # Finds PEC terms, [cm3/s]
                    (
                        cpls[upr][lwr]['EXCIT'][nd,nt],
                        cpls[upr][lwr]['RECOM'][nd,nt],
                        cpls[upr][lwr]['IONIZ'][nd,nt],
                        ) = crm.PECs(
                        # States for transition of interest
                        state = upr,
                        final = lwr,
                        # Data
                        data = data,
                        pops = pops,
                        nsim = nsim,   
                        # Plasma
                        ne_cm3 = ne_array[nd],
                        nz_cm3 = fz[0,:],
                        ZZ = ZZ,
                        )
                    ntrans += 1

    # Gets configuration names
    names = {}
    nn = 1
    for lvl in lvls:
        S2_1 = '2' ##################### Need to fix ####################

        names[lvl] = [
            str(nn),                                    # Index
            data[lvl]['sname'].decode("utf-8"),         # Configuration
            (
                str(nn)
                +'('+S2_1+')'
                +str(data[lvl]['P'])
                +'( {:1.1F})'.format(data[lvl]['2J']/2)
                )                                       # Quantum numbers
            ]

        # Increment
        nn += 1

    # Write adf15-format PEC file
    _write(
        path_fac=path_fac,
        cpls=cpls,
        ntrans=ntrans,
        names=names,
        ion=ion,
        ZZ=ZZ,
        Te_array=Te_array,
        ne_array=ne_array,
        )

###########################################
#
#           Utilities
#
###########################################

def _write(
    path_fac=None,
    cpls=None,
    ntrans=None,
    names=None,
    ion=None,
    ZZ=None,
    Te_array=None,
    ne_array=None,
    ):

    # Initialize file
    f = open(
        path_fac+'/pec#'+ion+str(ZZ)+'_cjp.dat',
        'w'
        )

    # File header
    f.write(
        "  "+str(ntrans)
        +"    /"+ion+" +"+str(ZZ)
        +" PHOTON EMISSIVITY COEFFICIENTS/\n"
        )
    f.write(
        "--------------------------------------------------------------------------------\n"
        )

    ##### ------ Data Section ------ #####

    # Initializes ISEL
    ISEL = 1

    # Initializes types
    types = ['IONIZ', 'EXCIT', 'RECOM']

    # Loop over types
    for ty in types:
        # Loop over upper states
        for upr in cpls.keys():
            # Loop over lower states
            for lwr in cpls[upr].keys():
                # Transition header
                f.write(
                    " {:1.6F}".format(cpls[upr][lwr]['lambda_A'])
                    +"   "+str(len(ne_array))
                    +"  "+str(len(Te_array))
                    +"    /TYPE = "+ty
                    +"   /ISEL = "+str(ISEL)
                    +"\n"
                    )

                # Writes data
                _write_block(
                    f=f,
                    data=cpls[upr][lwr][ty],
                    Te_array=Te_array,
                    ne_array=ne_array
                    )

                # Break
                f.write(
                    "--------------------------------------------------------------------------------\n"
                    )

                # Increments ISEL
                ISEL += 1

    ##### ------ Documentation Section ------ #####

    # Disclaimer
    f.write("C\n")
    f.write(
        "C PHOTON EMISSIVITY COEFFICIENT\n"
        )
    f.write(
        "C FROM FAC - USE WITH CAUTION\n"
        )
    f.write(
        "C NOT AN OFFICIAL ADAS PRODUCT!\n"
        )
    f.write("C\n")

    # Transitions header
    f.write(
        "C  ISEL  WAVELENGTH      TRANSITION               TYPE\n"
        )
    f.write(
        "C  ----  ----------  ---------------------------  -----\n"
        )

    ISEL = 1
    # Loop over types
    for ty in types:
        # Loop over upper states
        for upr in cpls.keys():
            # Loop over lower states
            for lwr in cpls[upr].keys():
                f.write(
                    "C   "+str(ISEL)+".  "
                    + "{:1.6F}".format(cpls[upr][lwr]['lambda_A'])
                    + names[upr][2].rjust(16, ' ')
                    + '-' + names[lwr][2].rjust(15, ' ')
                    + ty.rjust(7, ' ')
                    + "\n"
                    )

                ISEL +=1

    # Defining LS level configurations
    f.write("C\n")
    f.write(
        "C--------------------------------------------------------------------------------\n"
        )
    f.write("C\n")
    f.write("C Level configurations\n")
    f.write("C\n")

    for ind in np.arange(len(names)):
        f.write(
            "C"+names[ind][0].rjust(6, ' ')
            +". " +names[ind][1]
            +"\n"
            )

    # Signature
    f.write("C\n")
    f.write(
        "C--------------------------------------------------------------------------------\n"
        )
    f.write("C\n")

    f.write(
        "C    NOTES:\t PROCESSED \n"
        )
    f.write(
        "C\t\t Created by Conor Perks\n"
        )
    f.write(
        "C\t\t cjperks@psfc.mit.edu\n"
        )
    f.write("C\n")
    f.write(
        "C--------------------------------------------------------------------------------"
        )

    # Ends writing
    f.close()

#######################################################
#
#                 Utilities
#
#######################################################

def _write_block(
    f=None,
    data=None,      # dim(nd,nt)
    Te_array=None,  # dim(nt,)
    ne_array=None,  # dim(nd,)
    ):

    # Replaces zeros
    data[data<=1e-40] = 1e-40

    # Writes the density array block
    for rr in np.arange(int(np.ceil(len(ne_array)/8))):
        row = ''

        for cc in np.arange(8):
            try:
                row += (
                    "{:1.2E}".format(
                        ne_array[8*rr + ii]
                        ).rjust(9, ' ')
                    )
            except:
                blah = 0

        f.write(row+"\n")

    # Writes the temperature array block
    for rr in np.arange(int(np.ceil(len(Te_array)/8))):
        row = ''

        for cc in np.arange(8):
            try:
                row += (
                    "{:1.2E}".format(
                        Te_array[8*rr + ii]
                        ).rjust(9, ' ')
                    )
            except:
                blah = 0

        f.write(row+"\n")

    # Writes the data array block
    for nd in np.arange(len(ne_array)):
        for rr in np.arange(int(np.ceil(len(Te_array)/8))):
        row = ''

            for cc in np.arange(8):
                try:
                    row += (
                        "{:1.2E}".format(
                            data[nd, 8*rr + ii]
                            ).rjust(9, ' ')
                        )
                except:
                    blah = 0

            f.write(row+"\n")