'''

For the 2-photon emission distribution data table stored in ChiantiPy,
they kind of arbitrarily stopped at Z_nuclear=30 when the published
data goes up to Z_nuclear=92. So here, we fill the missing data

cjperks
Feb 19, 2024


'''

# Modules
from pfac import crm
import os
from scipy.interpolate import interp1d

# Paths
path2chianti = os.path.join(
    '/home/cjperks',
    'CHIANTI_10.1_database/continuum'
    )
'''
# Obtains original data files
f_ch_h = open(
    path2chianti+'/hseq_2photon.dat',
    'r'
    )
f_ch_he = open(
    path2chianti+'/heseq_2photon.dat',
    'r'
    )
'''
# Creates new data files
f_new_h = open(
    path2chianti+'/hseq_2photon_new.dat',
    'w'
    )
f_new_he = open(
    path2chianti+'/heseq_2photon_new.dat',
    'w'
    )

# Populates new data files with old data
for line in open(
    path2chianti+'/hseq_2photon.dat',
    'r'
    ):
    f_new_h.write(line)
for line in open(
    path2chianti+'/heseq_2photon.dat',
    'r'
    ):
    f_new_he.write(line)

######## Writes new data

import ChiantiPy.tools.io as io

dist_H = io.twophotonHRead() # []
dist_He = io.twophotonHeRead() # []

# Data on H-like distribution from
# S.P. Goldman & G.W.F. Drake, Phys. Rev. A, vol 24, num 1 (1981)
ys = np.r_[
    0.0,
    0.06250,
    0.12500,
    0.18750,
    0.25000,
    0.31250,
    0.37500,
    0.43750,
    0.50000
    ] # dim(ny,)
Zs = np.r_[
    1,
    20,	
    40,	
    60,	
    80,	
    92
    ] # dim(nZ,)

data = np.r_[
    0,0,0,0,0,0,
    2.03239,	1.94467,	1.71463,	1.41842,	1.12509,	0.96367,
    3.15792,	3.08515,	2.87808,	2.56769,	2.19190,	1.94784,
    3.84452,	3.78628,	3.61538,	3.34325,	2.98340,	2.72839,
    4.28435,	4.23831,	4.10084,	3.87404,	3.55733,	3.31977,
    4.56958,	4.53308,	4.42273,	4.23602,	3.96468,	3.75229,
    4.74855,	4.71885,	4.62821,	4.47188,	4.23760,	4.04805,
    4.84732,	4.82168,	4.74297,	4.60539,	4.39488,	4.22070,
    4.87892,	4.85464,	4.77990,	4.64866,	4.44629,	4.27750,
    ].reshape((len(ys), len(Zs))) # dim(ny, nZ)

y_mesh = np.concatenate(
    (
        ys, (1-np.flip(ys))[1:]
        )
    ) # dim(2ny-1,)

data_mesh = np.concatenate(
    (
        data,
        np.flip(data,axis=0)[1:, :]
        ),
    axis = 0
    ) # dim(2ny-1, nZ)

for zz in np.arange(31,92+1):
    if zz in Zs:
        ind = np.where(zz == Zs)[0][0]

        data_tmp = data_mesh[:,ind] # dim(2ny-1,)

    else:
        # Prepares interpolation
        # Assumes Z is high enough that He-like ~ H-like
        if zz <40:
            data_l_h = dist_H['psi0'][29]
            data_l_he = dist_He['psi0'][29]

            data_u = data_mesh[:,2]

            Z_l = 30
            Z_u = 40
        else:
            ind_l = np.where(Zs < zz)[0][-1]
            ind_u = np.where(Zs > zz)[0][0] 

            data_l_h = data_mesh[:,ind_l]
            data_l_he = data_mesh[:,ind_l]
            data_u = data_mesh[:,ind_u]

            Z_l = Zs[ind_l]
            Z_u = Zs[ind_u]

        # Interpolates between Z
        data_tmp = (
            (
                (data_u-data_l_h)
                /(Z_u-Z_l)
                )*(zz - Z_l)
            + data_l_h
            ) # dim(2ny-1,)

    # Normalization factor to area=2
    asum = np.trapz(data_tmp, y_mesh)/2

    # Prepares H-like text
    text_H = ''
    text_H += (
        "{:1.3F}".format(asum).rjust(10, ' ')
        )
    for ii in np.arange(len(y_mesh)):
        text_H += (
            "{:1.3F}".format(data_tmp[ii]).rjust(10, ' ')
            )

    f_new_h.write(
        str(zz).rjust(3, ' ')
        +" {:1.3e}".format(crm.TwoPhoton(zz,0))
        + text_H
        + '\n'
        )

    # Prepares He-like text
    data_tmp /= asum

    data_tmp_he = interp1d(
        y_mesh,
        data_tmp
        )(
            np.r_[[float(xx) for xx in dist_He['y0']]]
            )

    text_He = ''
    for ii in np.arange(len(dist_He['y0'])):
        text_He += (
            "{:1.3F}".format(data_tmp_he[ii]).rjust(10,' ')
            )

    f_new_he.write(
        str(zz).rjust(3, ' ')
        +" {:1.3e}".format(crm.TwoPhoton(zz,1))
        + text_He
        + '\n'
        )

# Closes files
f_new_h.close()
f_new_he.close()