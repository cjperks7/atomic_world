'''

Script to compare Mewe H-like cross-sections against
Mewe

cjperks

Jan 22, 223

'''

# Modules
from pfac import rfac
from atomic_world.calc_Mewe import EIE_fit as EIE
import os
import scipy.constants as cnt
from colradpy.convolve_EEDF import convolve_EEDF as cE

Te = np.logspace(np.log10(8.62e0), np.log10(8.62e4), 41) # [eV]

# Loads Mewe EIE data
Mewe = EIE.get_EIE(
    sp='Ar',
    Te = Te,
    nele = 1,
    )

# Loads FAC EIE data
paths = os.path.join(
    '/home/cjperks',
    '2008_SPARCxray/FAC/ATOM'
    )
vers = ['Ar', 'Ar_v2', 'Ar_v5', 'Ar_v6', 'Ar_v7']
ce = {}

fig, ax = plt.subplots(1,3)
colors = plt.cm.rainbow(np.linspace(0,1,len(vers)+1))

ax[0].plot(
    Mewe['Egrid'],
    Mewe['XS'][0,:],
    'r-',
    label = 'Mewe'
    )

ax[1].plot(
    Te,
    Mewe['PEC'][0,:]/1.06,
    'r-'
    )

ax[2].plot(
    Mewe['Egrid'],
    Mewe['XS'][0,:]/ (
        np.pi *cnt.physical_constants['Bohr radius'][0]**2
        * 13.6 /Mewe['Egrid']
        /(1 +2*0.5)
        *1e4
        ),
    'r-'
    )


for ii, ver in enumerate(vers):
    fil = os.path.join(
        paths,
        ver,
        'Ar01a.ce'
        )

    ce[ver] = {}
    ce[ver]['data'] = rfac.read_ce(fil)


    xx = (
        ce[ver]['data'][1][0]['born'][0,1]
        /(ce[ver]['data'][1][0]['Delta E'][0]*1e3)
        )
    ce[ver]['omega'] = (
        ce[ver]['data'][1][0]['bethe'][0]*np.log(xx) 
        + ce[ver]['data'][1][0]['born'][0,0]
        )
    ce[ver]['Eomega'] = ce[ver]['data'][1][0]['born'][0,1]

    #ce[ver]['XSomega'] = (
    #    np.pi * ce[ver]['omega']
    #    /(
    #        cnt.e * ce[ver]['data'][1][0]['born'][0,1]/cnt.hbar/cnt.c
    #        )**2
    #    /(1+ 2*0.5)
    #    ) *1e4
    ce[ver]['XSomega'] = (
        ce[ver]['omega']
        *np.pi *cnt.physical_constants['Bohr radius'][0]**2
        * 13.6 /ce[ver]['data'][1][0]['born'][0,1]
    /(1 +2*0.5)
    ) *1e4

    Etot = np.append(
        ce[ver]['data'][1][0]['EGRID'] + ce[ver]['data'][1][0]['TE0'],
        ce[ver]['Eomega']
        )

    XStot = np.append(
        ce[ver]['data'][1][0]['crosssection'][0,:]*1e-20,
        ce[ver]['XSomega']
        )

    ce[ver]['ratec'] = cE(
        EEDF = 'Maxwellian',
        Te=Te,
        XS = XStot,
        engyXS = Etot,
        m = 0,
        dE = None,
        )

    ax[2].plot(
        Etot,
        XStot /(
            np.pi *cnt.physical_constants['Bohr radius'][0]**2
            * 13.6 /Etot
            /(1 +2*0.5)
            *1e4
            ),
        '*',
        color = colors[ii]
        )

    ax[1].plot(
        Te,
        ce[ver]['ratec'],
        '*-',
        color = colors[ii]
        )


    ax[0].plot(
        ce[ver]['data'][1][0]['EGRID'] + ce[ver]['data'][1][0]['TE0'],
        ce[ver]['data'][1][0]['crosssection'][0,:]*1e-20,
        '*',
        label= ver,
        color = colors[ii]
        )

    ax[0].plot(
        ce[ver]['Eomega'],
        ce[ver]['XSomega'],
        'x',
        color = colors[ii]
        )



ax[0].set_yscale('log')
ax[0].set_xscale('log')
ax[0].set_xlabel('E_inc [eV]')
#ax.set_ylabel(r'$\Omega$')
ax[0].set_ylabel('XS [cm2]')
ax[0].grid('on')    
leg = ax[0].legend()
leg.set_draggable('on')

#ax[1].set_yscale('log')
#ax[1].set_xscale('log')
ax[1].set_xlabel('Te [eV]')
#ax.set_ylabel(r'$\Omega$')
ax[1].set_ylabel('ratec [cm3/s]')
ax[1].grid('on')    

ax[2].set_yscale('log')
ax[2].set_xscale('log')
ax[2].set_xlabel('E_inc [eV]')
ax[2].set_ylabel(r'$\Omega$')
ax[2].grid('on')    