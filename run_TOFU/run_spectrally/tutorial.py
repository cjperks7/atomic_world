'''

Tutorial script for spectrally

cjperks
July 16, 2024

'''

# Modules
import spectrally as sp


# ------------------------
# instanciate a Collection

coll = sp.Collection()

# -----------------
# create wavelength vector and gaussian
# -----------------

# ------------
# lamb

nlamb = 300
lamb = np.linspace(3.94e-10, 4.e-10, nlamb)

# ------------
# time vector

nt = 25
t = np.linspace(0, 10, nt)

# -----------
# spectra

lamb0 = 3.96e-10 + 0.001e-10 * np.cos(t[:, None])
lamb1 = 3.98e-10 + 0.001e-10 * np.sin(t[:, None])


# linear
linear = (1 + np.cos(t)[:, None]) * 10 * np.ones((1, nlamb))

gauss0 = (1 + np.cos(t)[:, None]) * 50 * np.exp(-(lamb[None, :]-lamb0)**2/0.002e-10**2)

gauss1= (1 + np.sin(t)[:, None]) * 80 * np.exp(-(lamb[None, :]-lamb1)**2/0.001e-10**2)


spectra = (
    linear
    + gauss0
    + gauss1
    + np.random.random((nt, nlamb)) * 5
)

# --------------------
# populate collection

# refs
coll.add_ref('nlamb', size=nlamb)
coll.add_ref('nt', size=nt)

# data - lamb
coll.add_data(
    'lamb',
    data=lamb,
    ref='nlamb',
    units='m',
)


# data - t
coll.add_data(
    't',
    data=t,
    ref='nt',
    units='s',
)


# data - spectrum
coll.add_data(
    'spectra',
    data=spectra,
    ref=('nt', 'nlamb'),
    units='counts',
)
dax = coll.plot_as_array('spectra', keyX='lamb', aspect='auto')


# ------------------
# add spectral model
# ------------------

# first mode : bck + 1 gauss + 1 lorentz
coll.add_spectral_model(
    key='sm0', 
    dmodel={
        'bck': 'exp_lamb', 
        'l0': {'type': 'gauss', 'lamb0': 3.96e-10}, 
        'l1': {'type': 'lorentz', 'lamb0': 3.7e-10},
    },
)
coll.show_details('sm0')


# first mode : bck + 1 gauss + 1 lorentz + 1 pvoigt
coll.add_spectral_model(
    key='sm1', 
    dmodel={
        'bck': 'linear', 
        'l0': {'type': 'gauss', 'lamb0': 3.96e-10}, 
        'l1': {'type': 'lorentz', 'lamb0': 3.7e-10},
        'l2': {'type': 'pvoigt', 'lamb0': 3.955e-10},
    },
    dconstraints={
        'g0': {
            'ref': 'l0_sigma',
            'l2_sigma': [0, 1, 0],
            'l1_gam': [0, 1, 0], 
            },
        'g2': {
            'ref': 'bck_a0',
            'bck_a1': [0, 0, 0],
        },
    },
)

# ------------------
# add spectral fit
# ------------------

coll.add_spectral_fit(
    key='sf0',
    key_model='sm0',
    key_data='spectra',
    key_lamb='lamb',
    dvalid={
        'nsigma': 3,
        'fraction': 0.41,
    },
)

dax = coll.plot_spectral_fit_input_validity()

coll.compute_spectral_fit()


# ------------------
# add atomic data
# ------------------

# Add line data from NIST
coll.add_spectral_lines_from_nist(
    lambmin=3.94e-10, 
    lambmax = 4e-10, 
    element='Ar'
    )

# Add line data from OpenADAS
coll.add_spectral_lines_from_openadas(
    lambmin=3.94e-10,
    lambmax=4e-10,
    element='Ar'
    )

# Add line data from external file
# NOTE: Example of formatting in:
# ~/miniconda3/envs/tofu3/lib/python3.9/site-packages/spectrally/tests/input/spectrallines.json
coll.add_spectral_lines_from_file()
