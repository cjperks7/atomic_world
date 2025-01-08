'''

Code to scan bend radius and thickness for the SPARC crystals

cjperks
Jan 8th, 2025

'''

# Modules
from pyTTE.pyTTE import TakagiTaupin, TTcrystal, TTscan, Quantity

# Crystal of interest
#cry = 'Kr'
#cry = 'Xe'
cry = 'LR'

# Default parameters
ddef = {
    'Kr':{
        'mat': 'Ge',
        'hkl': [2,4,2],
        'EE': 13.1, # [keV]
        'arc': np.linspace(-15,10,200),
        },
    'Xe':{
        'mat': 'AlphaQuartz',
        'hkl': [1,0,1],
        'EE': 4.56, # [keV]
        'arc': np.linspace(-30,20,200),
        },
    'LR':{
        'mat': 'Ge',
        'hkl': [2,0,2],
        'EE': 8.0, # [keV]
        'arc': np.linspace(-30,20,200),
        },
    }

# Scoping parameters
RR_sag = np.linspace(0.1, 1, 10) # [m]
R0_sag = 0.5 # [m]
if cry == 'LR':
    t0 = 0.1 # [mm]
    #tt = np.linspace(0.02, 0.2, 10)
    tt = np.linspace(0.1, 0.7, 10)
else:
    t0 = 0.5 # [mm]
    tt = np.linspace(0.1, 1, 10) # [mm]
RR_mer = np.linspace(-1, -0.1, 10) # [m]
R0_mer = -0.5
#RR_mer = np.linspace(0.1, 1, 10)
#R0_mer = 0.5


# Defines crystal
ttx_flat = TTcrystal(
    crystal = ddef[cry]['mat'],     # Material
    hkl = ddef[cry]['hkl'],         # Miller index
    thickness = Quantity(t0, 'mm'),  # thickness
    debye_waller = 1,               # ???
    Rx = 'inf',                     # Bend radius in meridional dir
    Ry = 'inf',                     # Bend radius in saggital dir
    )

# Angular scan for rocking curve
tts_sig = TTscan(
    constant=Quantity(ddef[cry]['EE'], 'keV'),
    scan = Quantity(ddef[cry]['arc'], 'arcsec'),
    polarization = 'sigma'
    )
tts_pi = TTscan(
    constant=Quantity(ddef[cry]['EE'], 'keV'),
    scan = Quantity(ddef[cry]['arc'], 'arcsec'),
    polarization = 'pi'
    )

# Simulates a flat crystal
tt_flat_s = TakagiTaupin(ttx_flat, tts_sig)
tt_flat_p = TakagiTaupin(ttx_flat, tts_pi)

vec_s, R_s, T_s = tt_flat_s.run()
vec_p, R_p, T_p = tt_flat_p.run()

vec_flat = vec_s.copy()
R_flat = 0.5*(R_s+R_p)

# Init figure
if cry == 'LR':
    fig, ax = plt.subplots(2,3)
else:
    fig, ax = plt.subplots(2,2)


ax[0,0].plot(
    vec_flat, 
    R_flat, 
    '-',
    color = 'r',
    label = 'flat'
    )
ax[1,0].plot(
    vec_flat - vec_flat[np.argmax(R_flat)],
    R_flat,
    '-',
    color = 'r',
    label = 'flat'
    )

for rad in RR_sag:
    # Defines crystal
    if cry == 'LR':
        ttx_bent = TTcrystal(
            crystal = ddef[cry]['mat'],     # Material
            hkl = ddef[cry]['hkl'],         # Miller index
            thickness = Quantity(t0, 'mm'),  # thickness
            debye_waller = 1,               # ???
            Rx = Quantity(R0_mer, 'm'),     # Bend radius in meridional dir
            Ry = Quantity(rad, 'm'),        # Bend radius in saggital dir
            )
    else:
        ttx_bent = TTcrystal(
            crystal = ddef[cry]['mat'],     # Material
            hkl = ddef[cry]['hkl'],         # Miller index
            thickness = Quantity(t0, 'mm'),  # thickness
            debye_waller = 1,               # ???
            Rx = 'inf',                     # Bend radius in meridional dir
            Ry = Quantity(rad, 'm'),        # Bend radius in saggital dir
            )

    # Simulates a flat crystal
    tt_bent_s = TakagiTaupin(ttx_bent, tts_sig)
    tt_bent_p = TakagiTaupin(ttx_bent, tts_pi)

    vec_s, R_s, T_s = tt_bent_s.run()
    vec_p, R_p, T_p = tt_bent_p.run()
    R_bent = 0.5*(R_s+R_p)

    ax[0,0].plot(
        vec_s,
        R_bent,
        '-',
        label = 'Ry=%0.1f m'%(rad)
        )
    ax[1,0].plot(
        vec_s - vec_s[np.argmax(R_bent)],
        R_bent,
        '-',
        label = 'Ry=%0.1f m'%(rad)
        )

ax[0,0].grid('on')
ax[0,0].set_xlabel(r'$\theta-\theta_B$ [arcsec]')
ax[0,0].set_ylabel('Reflected power [frac]')
if cry == 'LR':
    ax[0,0].set_title(r'Bend radius (saggital) scan (t=%0.1f mm, R$_{mer}$ = %0.1f m)'%(t0, R0_mer))
else:
    ax[0,0].set_title('Bend radius (saggital) scan (t=%0.1f mm)'%(t0))
leg = ax[0,0].legend(labelcolor='linecolor')
leg.set_draggable('on')

ax[1,0].grid('on')
ax[1,0].set_xlabel(r'$\theta-\theta_{peak}$ [arcsec]')
ax[1,0].set_ylabel('Reflected power [frac]')




ax[0,1].plot(
    vec_flat, 
    R_flat, 
    '-',
    color = 'r',
    label = 'flat'
    )
ax[1,1].plot(
    vec_flat - vec_flat[np.argmax(R_flat)],
    R_flat,
    '-',
    color = 'r',
    label = 'flat'
    )

for thi in tt:
    # Defines crystal
    if cry == 'LR':
        ttx_bent = TTcrystal(
            crystal = ddef[cry]['mat'],     # Material
            hkl = ddef[cry]['hkl'],         # Miller index
            thickness = Quantity(thi, 'mm'),  # thickness
            debye_waller = 1,               # ???
            Rx = Quantity(R0_mer, 'm'),     # Bend radius in meridional dir
            Ry = Quantity(R0_sag, 'm'),        # Bend radius in saggital dir
            )
    else:
        ttx_bent = TTcrystal(
            crystal = ddef[cry]['mat'],     # Material
            hkl = ddef[cry]['hkl'],         # Miller index
            thickness = Quantity(thi, 'mm'),  # thickness
            debye_waller = 1,               # ???
            Rx = 'inf',                     # Bend radius in meridional dir
            Ry = Quantity(R0_sag, 'm'),        # Bend radius in saggital dir
            )

    # Simulates a flat crystal
    tt_bent_s = TakagiTaupin(ttx_bent, tts_sig)
    tt_bent_p = TakagiTaupin(ttx_bent, tts_pi)

    vec_s, R_s, T_s = tt_bent_s.run()
    vec_p, R_p, T_p = tt_bent_p.run()
    R_bent = 0.5*(R_s+R_p)

    ax[0,1].plot(
        vec_s,
        R_bent,
        '-',
        label = 't=%0.1f mm'%(thi)
        )
    ax[1,1].plot(
        vec_s - vec_s[np.argmax(R_bent)],
        R_bent,
        '-',
        label = 't=%0.1f mm'%(thi)
        )

ax[0,1].grid('on')
ax[0,1].set_xlabel(r'$\theta-\theta_B$ [arcsec]')
ax[0,1].set_ylabel('Reflected power [frac]')
if cry == 'LR':
    ax[0,1].set_title(r'Thickness scan (R$_{sag}$=%0.1f m, R$_{mer}$=%0.1f m)'%(R0_sag, R0_mer))
else:
    ax[0,1].set_title('Thickness scan (R=%0.1f m)'%(R0_sag))
leg = ax[0,1].legend(labelcolor='linecolor')
leg.set_draggable('on')

ax[1,1].grid('on')
ax[1,1].set_xlabel(r'$\theta-\theta_{peak}$ [arcsec]')
ax[1,1].set_ylabel('Reflected power [frac]')




if cry == 'LR':
    ax[0,2].plot(
        vec_flat, 
        R_flat, 
        '-',
        color = 'r',
        label = 'flat'
        )
    ax[1,2].plot(
        vec_flat - vec_flat[np.argmax(R_flat)],
        R_flat,
        '-',
        color = 'r',
        label = 'flat'
        )

    for rad in RR_mer:
        # Defines crystal
        ttx_bent = TTcrystal(
            crystal = ddef[cry]['mat'],     # Material
            hkl = ddef[cry]['hkl'],         # Miller index
            thickness = Quantity(t0, 'mm'),  # thickness
            debye_waller = 1,               # ???
            Rx = Quantity(rad, 'm'),     # Bend radius in meridional dir
            Ry = Quantity(R0_sag, 'm'),        # Bend radius in saggital dir
            )

        # Simulates a flat crystal
        tt_bent_s = TakagiTaupin(ttx_bent, tts_sig)
        tt_bent_p = TakagiTaupin(ttx_bent, tts_pi)

        vec_s, R_s, T_s = tt_bent_s.run()
        vec_p, R_p, T_p = tt_bent_p.run()
        R_bent = 0.5*(R_s+R_p)

        ax[0,2].plot(
            vec_s,
            R_bent,
            '-',
            label = 'Ry=%0.1f m'%(rad)
            )
        ax[1,2].plot(
            vec_s - vec_s[np.argmax(R_bent)],
            R_bent,
            '-',
            label = 'Ry=%0.1f m'%(rad)
            )

    ax[0,2].grid('on')
    ax[0,2].set_xlabel(r'$\theta-\theta_B$ [arcsec]')
    ax[0,2].set_ylabel('Reflected power [frac]')
    ax[0,2].set_title(r'Bend radius (meridional) scan (t=%0.1f mm, R$_{sag}$ = %0.1f m)'%(t0, R0_sag))
    leg = ax[0,2].legend(labelcolor='linecolor')
    leg.set_draggable('on')

    ax[1,2].grid('on')
    ax[1,2].set_xlabel(r'$\theta-\theta_{peak}$ [arcsec]')
    ax[1,2].set_ylabel('Reflected power [frac]')



fig.suptitle(
    ddef[cry]['mat']
    +'(%i%i%i)'%(ddef[cry]['hkl'][0], ddef[cry]['hkl'][1], ddef[cry]['hkl'][2])
    +r'; $E_x$ = %0.2f keV'%(ddef[cry]['EE'])
    )