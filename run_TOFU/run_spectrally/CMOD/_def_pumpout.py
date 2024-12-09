'''

Stores default settings for analyzing VUV spectra
for the pumpout shots

cjperks
Aug 30, 2024

'''

# Modules
import numpy as np
import scipy.constants as cnt

__all__ = [
    '_def_slicing',
    '_def_dsettings'
    ]

###############################
#
#           Main
#
###############################

# Lambda and time slicing for VUV spectra
def _def_slicing(
    shot =None,
    system = None,
    ):

    ds = {
        1140221013: {
            'XEUS' :{
                'lamb_AA': [23, 29], # [AA]
                'time_s': [0.84, 1.44], # [s]
                },
            'LoWEUS' :{
                'lamb_AA': [110, 230], # [AA]
                'time_s': [0.84, 1.44], # [s]
                }
            },
        1120913016: {
            'XEUS' :{
                'lamb_AA': [23, 29], # [AA]
                'time_s': [0.84, 1.44], # [s]
                },
            'LoWEUS' :{
                'lamb_AA': [110, 230], # [AA]
                'time_s': [0.84, 1.44], # [s]
                }
            },
        }

    # Output
    return ds[shot][system]

# Spectral fitting settings
def _def_dsettings(
    shot = None,
    system = None,
    ):

    dset = {
        1140221013: {
            'XEUS': {
                'dmodel': {
                    'Ar':{ # Ar lines of interest
                        'mz': np.ones(11)*39.948*cnt.m_u, # [kg]
                        'vlamb': [
                            23.73e-10, # Ar XVI
                            24.97e-10, # Ar XV
                            25.09e-10, # Ar XVI
                            25.21e-10, # Ar XVI
                            25.75e-10, # Ar XVI
                            25.92e-10, # Ar XVI
                            26.16e-10, # Ar XV
                            27.28e-10, # Ar XV
                            27.67e-10, # Ar XIV
                            27.87e-10, # Ar XIV
                            28.58e-10, # Ar XV
                            ],
                        },
                    'Misc':{ # Unidentified lines
                        'mz': np.ones(18)*95.95*cnt.m_u, # [kg]
                        'vlamb': [
                            23.05e-10,
                            23.18e-10,
                            23.48e-10,
                            23.88e-10,
                            24.08e-10,
                            24.25e-10,
                            24.83e-10,
                            25.42e-10,
                            26.27e-10,
                            26.55e-10,
                            26.76e-10,
                            26.87e-10,
                            27.43e-10,
                            28.13e-10,
                            28.33e-10,
                            28.95e-10,
                            26.02e-10,
                            28.03e-10,
                            ]
                        },
                    'bck': 'poly'
                    },
                'dconstraints': {
                    'g0': {
                        'ref': 'bck_a0',
                        'bck_a1': [0,0,0], # A*bck_a0**0 + B*bck_a0**1 + C*bck_a0**2
                        'bck_a2': [0,0,0], # A*bck_a0**0 + B*bck_a0**1 + C*bck_a0**2
                        }, # a0*lambda**0 + a1*lambda**1 +a2*lambda**2
                    },
                'dvalid': {
                    'nsigma': 0.01, # minimum signal-to-noise ratio
                    'fraction': 0.1, # fraction of valid pixels in spectrum for the whole thing to be valid
                    },
                'dbnds': {
                    'Ar':{
                        'dlamb': 3.9e-12, # [m]
                        'vc_bnd': 0.02,
                        'sig_bnd': 1.1
                        },
                    'Misc':{
                        'dlamb': 3.9e-12, # [m]
                        'vc_bnd': 0.02,
                        'sig_bnd': 5
                        }
                    }
                },
            'LoWEUS': {
                'dmodel': {
                    'Ar': {
                        'mz': [
                            39.948*cnt.m_u
                            ], # [kg]
                        'vlamb': [
                            224.9e-10, # Ar XV
                            ],
                        },
                    'Mo': {
                        'mz': [
                            95.95*cnt.m_u,
                            95.95*cnt.m_u,
                            ], # [kg]
                        'vlamb': [
                            118.7e-10, # Mo XXXI
                            130.7e-10, # Mo XXXII
                            ],
                        },
                    'bck': 'poly'
                    #'bck': 'exp_lamb'
                    },
                'dconstraints': {
                    'g0': {
                        'ref': 'bck_a0',
                        'bck_a1': [0,0,0], # A*bck_a0**0 + B*bck_a0**1 + C*bck_a0**2
                        'bck_a2': [0,0,0], # A*bck_a0**0 + B*bck_a0**1 + C*bck_a0**2
                        }, # a0*lambda**0 + a1*lambda**1 +a2*lambda**2
                    },
                #'dconstraints': None,
                'dvalid': {
                    'nsigma': 0.01, # minimum signal-to-noise ratio
                    'fraction': 0.1, # fraction of valid pixels in spectrum for the whole thing to be valid
                    },
                'dbnds': {
                    'Ar':{
                        'dlamb': 3.9e-11, # [m]
                        'vc_bnd': 0.02,
                        'sig_bnd': 1.5
                        },
                    'Mo':{
                        'dlamb': 3.9e-11, # [m]
                        'vc_bnd': 0.02,
                        'sig_bnd': 1.5
                        },
                    }
                },   
            },
        }

    # Output
    return dset[shot][system]