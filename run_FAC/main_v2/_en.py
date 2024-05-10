'''

_en.py is a script meant to facilitate calculating
atomic energy levels using FAC

cjperks
May 9, 2024

'''

# Modules
import numpy as np
import sys, os
import copy
import time

__all__ = [
    '_states',
    '_run_en'
    ]

######################################################
#
#                   Main
#
######################################################

# Generates configuration groups
def _states(
    name=None,
    nele=None,
    settings=None,
    fac = None,
    ):

    print('Generating states of interest')

    groups = {
        'grd': [],
        'exc': [],
        'ion': [],
        }


     # Loop over group types
    for typ in groups.keys():
        print(typ)

        # If considering a ground state
        if typ == 'grd':
            for ii in np.arange(len(settings['grd']['1'])):
                # Builds ground state
                groups = _build_state(
                    ind = ii,
                    label = ii,
                    config = settings['grd'],
                    typ = typ,
                    groups = groups,
                    fac = fac,
                    restrict = settings['restrictions'],
                    )

        # If considering an excited state
        elif typ == 'exc':
            grd_ind = settings['single']['grd_ind']

            # ---- Single excited electron --- #
            print('single')

            # Loop over electrons to excite
            for ii in np.arange(len(settings['single']['n_raise'])):
                # Electron to delete
                nd = settings['single']['n_raise'][ii]

                # Loop levels to excite to
                for nn in range(
                    settings['single']['n_min'][ii], 
                    settings['single']['n_max'][ii]+1
                    ):
                    tmp = copy.deepcopy(settings['grd'])

                    # Deletes electron to excite
                    tmp[str(nd)][grd_ind] -= 1

                    # Promotes electron
                    if str(nn) not in tmp.keys():
                        tmp[str(nn)] = [0]*len(tmp['1'])
                    tmp[str(nn)][grd_ind] += 1

                    # Builds ground state
                    groups = _build_state(
                        ind = grd_ind,
                        label = [nd,nn],
                        config = tmp,
                        typ = typ,
                        groups = groups,
                        fac = fac,
                        long = True,
                        restrict = settings['restrictions'],
                        )

            # ---- Double excited electrons --- #
            print('double')

            # Loop over electrons to excite
            for ii in np.arange(len(settings['double']['n_raise'])):
                # Electron to delete
                nd1 = settings['double']['n_raise'][ii][0]
                nd2 = settings['double']['n_raise'][ii][1]

                # Loop over levels to excite to
                for nn1 in range(
                    settings['double']['n_min'][ii][0], 
                    settings['double']['n_max'][ii][0]+1
                    ):
                    for nn2 in range(
                        settings['double']['n_min'][ii][1], 
                        settings['double']['n_max'][ii][1]+1
                        ):
                        tmp = copy.deepcopy(settings['grd'])

                        # Avoid double counting
                        if nn2 < nn1:
                            continue

                        # Deletes electron to excite
                        tmp[str(nd1)][grd_ind] -= 1
                        tmp[str(nd2)][grd_ind] -= 1

                        # Promotes electron
                        if str(nn1) not in tmp.keys():
                            tmp[str(nn1)] = [0]*len(tmp['1'])
                        if str(nn2) not in tmp.keys():
                            tmp[str(nn2)] = [0]*len(tmp['1'])
                        tmp[str(nn1)][grd_ind] += 1
                        tmp[str(nn2)][grd_ind] += 1

                        # Builds ground state
                        groups = _build_state(
                            ind = grd_ind,
                            label = [nd1, nd2,nn1,nn2],
                            config = tmp,
                            typ = typ,
                            groups = groups,
                            fac = fac,
                            double = True,
                            long = True,
                            restrict = settings['restrictions'],
                            )

        # If considering an ionized state
        elif typ == 'ion':
            # Loop over electrons to ionize
            for ii in np.arange(len(settings['ion']['n_remove'])):
                tmp = copy.deepcopy(settings['grd'])

                # Electron to delete
                nd = settings['ion']['n_remove'][ii]

                # Deletes electron to excite
                tmp[str(nd)][grd_ind] -= 1

                # Builds ground state
                groups = _build_state(
                    ind = grd_ind,
                    label = ii,
                    config = tmp,
                    typ = typ,
                    groups = groups,
                    fac = fac,
                    restrict = settings['restrictions'],
                    )

    print(groups)
    # Output
    return groups


# Generates energy levels
def _run_en(
    name = None,
    nele = None,
    groups = None,
    settings = None,
    fac = None,
    ):

    # Track time spent
    start = time.time()
    start0 = start

    print('EN: Calc energy levels')

    # Loop over groups
    for typ in groups.keys():
        # If considering a ground states
        if typ == 'grd':
            print('Structure: ground complex')
            if (nele > 1):
                fac.ConfigEnergy(0)
            fac.OptimizeRadial(groups['grd'][0])
            if (nele > 1):
                fac.ConfigEnergy(1)

            structs = [xx for xx in groups[typ]]
            structs.append(groups['exc'][0])
            print(structs)

            # Calculates atomic structure with first excited state
            fac.Structure(name+'b.en', structs)

            # Precalculates angular coefficients with first excited state
            fac.PrepAngular(structs)

        # If considering an ionized state
        elif typ == 'ion':
            print('Structure: ionized complex')

            for ii,ion in enumerate(groups[typ]):
                print([ion])
                # Calculates atmoic structure
                fac.Structure(name+'b.en', [ion])

                # Precalculates angular coefficients
                fac.PrepAngular([ion])
                for jj in range(ii+1,len(groups[typ])):
                    fac.PrepAngular([ion], [groups[typ][jj]])

                fac.PrepAngular([ion], [groups['grd'][0]])

        # If considering an excited state
        elif typ == 'exc':
            print('Structure: excited complexes')

            # Loop over excited states excluding the first
            for st in groups[typ][1:]:
                print([st])
                # Calculates atomic structure
                fac.Structure(name+'b.en', [st])

                # Precalculates angular coefficients
                fac.PrepAngular([groups['ion'][0]], [st])
                fac.PrepAngular([groups['grd'][0]], [st])


    # Constructs spectator electron for recombining states
    fac.SetRecPWOptions(settings['spectator']['rec_pw_max'])
    fac.SetRecSpectator(
        settings['spectator']['n_spec_min'],
        settings['spectator']['n_frozen']
        )

    st = groups['ion'][settings['spectator']['ion_ind']]
    for nn in range(
        settings['spectator']['n_spec_min'], 
        settings['spectator']['n_spec_max']+1
        ):
        print(([st], nn))
        fac.RecStates(name+'b.en',
            [st],
            nn
            )
    # Highly excited state
    if settings['spectator']['n_ext'] is not None:
        print(([st], settings['spectator']['n_ext']))
        fac.RecStates(name+'b.en',
            [st],
            settings['spectator']['n_ext']
            )


    # Print the time this took
    stop = time.time()
    s='Done %10.3E s'%(stop-start)
    print(s)


######################################################
#
#                   Utilities
#
######################################################

# Builds dictionary to state string
def _build_state(
    ind = None,
    label = None,
    config = None,
    typ = None,
    groups = None,
    fac = None,
    double = False,
    long = False,
    restrict = None,
    ):

    # Init structure string
    struct = ''

    for nn in config.keys():
        if config[nn][ind] == 0:
            continue
        else:
            struct += '%i*%i '%(int(nn),config[nn][ind])
    struct = struct[:-1]

    # Includes restrictions
    struct += _add_restrict(config=config, restrict=restrict)

    # Defines ground state
    if long:
        if double:
            groups[typ].append(typ+'.%i.%i->%i.%i'%(
                label[0], label[1], label[2], label[3]
                ))
        else:
            groups[typ].append(typ+'.%i->%i'%(
                label[0], label[1]
                ))
    else:
        groups[typ].append(typ+'.%i'%(label))

    fac.Config(groups[typ][-1], struct)
    print(struct)

    # Output
    return groups


# Add level restrictions
def _add_restrict(
    config = None,
    restrict = None,
    ):

    # Init
    out = ''

    # Orbitals
    orb = {
        '0': 's',
        '1': 'p',
        '2': 'd',
        '3': 'f',
        '4': 'g',
        '5': 'h',
        '6': 'i',
        '7': 'k',
        '8': 'l',
        '9': 'm',
        '10': 'n'
        }

    # If don't add any
    if restrict is None:
        return out

    elif 'l_max' in restrict.keys():

        # Loop over principal quantum number
        for nn in config.keys():
            if int(nn) > restrict['l_max']:
                for ii in range(restrict['l_max'], int(nn)):
                    out += (
                        ';'
                        + str(nn)
                        + orb[str(ii)]
                        + '<1'
                        )

    # Output
    return out
