'''

Trying out the pyTTE package for rocking curves with bent crystals

cjperks
Jan 7th, 2025

'''

# Modules
from pyTTE.pyTTE import TakagiTaupin, TTcrystal, TTscan, Quantity

# Defines crystal
ttx_Kr = TTcrystal(
    crystal = 'Ge',                 # Material
    hkl = [2,4,2],                  # Miller index
    thickness = Quantity(0.1, 'mm'),  # thickness
    debye_waller = 1,               # ???
    Ry = 'inf',#Quantity(0.6, 'm'),        # Bend radius in saggital dir
    Rx = 'inf',                     # Bend radius in meridional dir
    fix_to_axes = 'shape',          # How crystal is bent: 'shape' --> glue to substrate, 'torque' --> edge pins
    )

ttx_Kr_flat = TTcrystal(
    crystal = 'Ge',                 # Material
    hkl = [2,4,2],                  # Miller index
    thickness = Quantity(1, 'mm'),  # thickness
    debye_waller = 1,               # ???
    Rx = 'inf',                     # Bend radius in one dir
    Ry = 'inf',                     # Bend radius in other dir
    )

ttx_Xe = TTcrystal(
    crystal = 'AlphaQuartz',                 # Material
    hkl = [1,0,1],                  # Miller index
    thickness = Quantity(1, 'mm'),  # thickness
    debye_waller = 1,               # ???
    Rx = Quantity(0.6, 'm'),        # Bend radius in one dir
    Ry = 'inf',                     # Bend radius in other dir
    )

ttx_Xe_flat = TTcrystal(
    crystal = 'AlphaQuartz',                 # Material
    hkl = [1,0,1],                  # Miller index
    thickness = Quantity(1, 'mm'),  # thickness
    debye_waller = 1,               # ???
    Rx = 'inf',                     # Bend radius in one dir
    Ry = 'inf',                     # Bend radius in other dir
    )

# Angular scan for rocking curve
tts_Kr_s = TTscan(
    constant=Quantity(13.1, 'keV'),
    scan = Quantity(np.linspace(-5,10,150), 'arcsec'),
    polarization = 'sigma'
    )
tt_Kr_flat_s = TakagiTaupin(ttx_Kr_flat, tts_Kr_s)
tt_Kr_s = TakagiTaupin(ttx_Kr, tts_Kr_s)

tts_Xe_s = TTscan(
    constant=Quantity(4.56, 'keV'),
    scan = Quantity(np.linspace(-10,30,150), 'arcsec'),
    polarization = 'sigma'
    )
tt_Xe_flat_s = TakagiTaupin(ttx_Xe_flat, tts_Xe_s)
tt_Xe_s = TakagiTaupin(ttx_Xe, tts_Xe_s)



ddata = {}
ddata['Kr'] = {}
ddata['Xe'] = {}

ddata['Kr']['flat'] = {}
ddata['Kr']['bent'] = {}

scan_vector, R, T = tt_Kr_flat_s.run()
ddata['Kr']['flat'] = {'angle': scan_vector, 'pwr': R}

scan_vector, R, T = tt_Kr_s.run()
ddata['Kr']['bent']= {'angle': scan_vector, 'pwr': R}




scan_vector, R, T = tt_Xe_flat_s.run()
ddata['Xe']['flat'] = {'angle': scan_vector, 'pwr': R}

scan_vector, R, T = tt_Xe_s.run()
ddata['Xe']['bent']= {'angle': scan_vector, 'pwr': R}



sym = 'Kr'
fig, ax = plt.subplots()

ax.plot(
    ddata[sym]['flat']['angle'],
    ddata[sym]['flat']['pwr'],
    color = 'r',
    label = 'flat'
    )

ax.plot(
    ddata[sym]['bent']['angle'],
    ddata[sym]['bent']['pwr'],
    color = 'b',
    label = 'bent'
    )

ax.grid('on')
leg = ax.legend(labelcolor='linecolor')
leg.set_draggable('on')

ax.set_xlabel(r'$\theta-\theta_B$ [arcsec]')
ax.set_ylabel('Reflected power')





ttx = TTcrystal(
    crystal='Si', 
    hkl=[6,6,0], 
    thickness=Quantity(150,'um'), 
    R = Quantity(100,'cm'), 
    asymmetry = Quantity(5, 'deg'), 
    nu = 0.22, 
    E = Quantity(150, 'GPa')
    )
tts = TTscan(
    scan = 200, 
    constant = Quantity(75, 'deg'), 
    polarization = 'pi'
    )

tt = TakagiTaupin(ttx,tts)
tt.run()
tt.plot()