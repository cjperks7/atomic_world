'''

Script to use ChiantiPy to calc lines in VUV region

cjperks
May 22, 2024

'''

# Modules
import ChiantiPy.core as ch

# Plasma
Te_K = 1e6 # [K]
ne_cm3 = 1e12 # [1/cm3]

sps = [
    #'C',
    #'O',
    #'N'
    'Fe'
    ]

css = [
    #[4,5],
    #[4,5,6,7]
    #[4,5,6]
    #[9,10,12,13]
    [14]
    ] # spectroscopic notation

ylim = [40,70] # [AA]

for ii, sp in enumerate(sps):
    for cs in css[ii]:
        ion = ch.ion(
            sp.lower()+'_%i'%(cs),
            temperature = Te_K,
            eDensity = ne_cm3,
            em = 1
            )

        ion.intensityPlot(
            wvlRange=ylim
            )