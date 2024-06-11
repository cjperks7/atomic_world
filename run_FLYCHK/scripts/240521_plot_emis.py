'''

Script to plot spectra caluclated using FLYCHK

cjperks
May 21, 2024

'''

# Modules
from atomic_world.run_FLYCHK.read import read_emis as re

#Te_eVs = np.linspace(5e3,9.5e3,10)
Te_eVs = [9.5e3]


for Te_eV in Te_eVs:
    # Loads spectra
    spec = re.read(
        #suffix = '_E4200-4800eV',
        suffix = '_E1000-6000eV',
        #suffix = '_E6000-9000eV',
        sp = 'W',
        Te_eV = Te_eV,
        ne_cm3 = 1e14,
        #Ti_eV = 2e4,
        nimp_cm3 = 1.5e-5*1e14,
        convolve = True
        )

    # Plots spectra
    re.plot_spec(
        out=spec,
        sp = 'W',
        x_units = 'eV',
        incl_cs = [12,13,14,15,16,17,18,19,20],
        xlim = np.r_[2.7,2.8], # [AA]
        #xlim = np.r_[4.13,12.3],
        tol = 1e-2,
        calib_vals = np.r_[2.7490, 2.7520, 2.72], # [AA]
        calib_wdth = 0.0005 # [AA]
        )





    import scipy.constants as cnt
    hc = cnt.h*cnt.c/cnt.e*1e10 # [eV*AA]
    sp = 'W'

    for cs in spec[sp].keys():
        with open(
            '/home/cjperks/test/W_cjperks_FLYCHK_spec_Te%0.1fkeV_%i_E1000-6000eV.txt'%(Te_eV/1e3,cs),
            'w') as file:
            # Zip the arrays and write each pair to the file
            for item1, item2 in zip(
                hc/spec[sp][cs]['spec_lambda_A'], 
                spec[sp][cs]['spec_ph/cm3/s/AA']
                ):
                file.write(f"{item1}\t{item2}\n")