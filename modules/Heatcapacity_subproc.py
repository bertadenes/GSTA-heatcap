#!/usr/bin/env python3
"""
This part does the calculation of heat capacity.
"""

from modules.GetHeatcapacity_util import *
from multiprocessing import Pool
import modules.g09BOMD_filter as g09f
from functools import partial
from contextlib import closing
import warnings
import numpy as np
from os import listdir
from os import getcwd

try:
    import readline

    rl = True
except ImportError:
    rl = False
    pass

warnings.filterwarnings('ignore')
try:
    import matplotlib.pyplot as plt
    plot = True
except ImportError:
    # print("\nFailed importing matplotlib \nNo plotting, sorry")
    plot = False

# This part is needed for tab autocomplete thing, when code is run interactively
ldir = listdir(getcwd())


def completer(text, state):
    options = [x for x in ldir if x.startswith(text)]
    try:
        return options[state]
    except IndexError:
        return None


if rl:
    readline.set_completer(completer)
    readline.parse_and_bind("tab: complete")
else:
    pass


def DOS_from_velocity(path_pos, path_neg, calcdos=True, plot=False, temp=None):
    """
    Function to compute density of state function via the Fourier transformation
    of the velocity autocorrelation function 
    :param path_pos: Path to positive BOMD output
    :param path_neg: Path to negative BOMD output
    :param calcdos: True by default
    :param plot: True if plotting of functions is required (False by default) 
    :param temp: Temperature of smoothing
    :return: Density of states function
    """
    J2Eh = 2.293710449e17  # Joule to Hartree conversion
    Eh2J = J2Eh ** -1  # Hartree to Joule conversion
    MwVelp = []
    MwVeln = []
    MwVel = []
    if temp is None:
        temp = g09f.getTemp(path_pos)
    AtMass = g09f.getAtMass(path_pos)

    # Extract mass weighted velocities from outputs
    for j in np.nditer(g09f.getVel(path_pos), flags=['external_loop'], order='F'):
        MwVelp.append(j)
    for j in np.nditer(g09f.getVel(path_neg), flags=['external_loop'], order='F'):
        MwVeln.append(j)

    # Concatenate positive and negative arrays
    j = 0
    for i in MwVeln:
        MwVel.append(np.concatenate((np.flipud(i), MwVelp[j][1:])))
        j += 1
    j = 0

    Na, Nt, Ns, dt = g09f.getCalcInfo(path_pos)
    t = np.empty(len(MwVel[0]))

    # get time array
    j = 0
    i = 0
    while i < len(t):
        # cartvel[0][i] = np.cos(0.1*i)# for test only
        t[i] = j
        j += dt * 1e-15
        i += 1

    # Smoothing velocities
    # smooth = partial(g09f.smoothFn_windowed, arg=t, Tf=temp, target=2, a=None, b=None)
    # with closing(Pool(processes=8)) as pool:
    #     smvel = pool.map(smooth, MwVel)
    #     pool.terminate()

    # Computing Density of state function
    if calcdos:
        acoriginal = []
        for j in MwVel:
            acoriginal.append(autocorr_manual(j))
        actotone = sumfn(acoriginal)
        dosoriginal = getDos(actotone)

        # ac = []
        # for i in smvel:
        #     ac.append(autocorr_manual(i))
        # acsmvel = sumfn(ac)
        # dossmvel = getDos(acsmvel)
    else:
        pass

    #### Plotting some of the functions if required ####
    if plot and calcdos:
        x1 = range(6000)
        plt.figure(1)
        plt.subplot(211)
        plt.plot(x1, dosoriginal)
        # plt.subplot(212)
        # plt.plot(x1, dossmvel)
        # plt.show()
    else:
        pass

    return dosoriginal  # , dossmvel


def DOS_from_velocity_single(path_pos, calcdos=True, plot=False, temp=None):
    """
    Function to compute density of state function via the Fourier transformation
    of the velocity autocorrelation function
    :param path_pos: Path to positive BOMD output
    :param calcdos: True by default
    :param plot: True if plotting of functions is required (False by default)
    :param temp: Temperature of smoothing
    :return: Density of states function
    """
    J2Eh = 2.293710449e17  # Joule to Hartree conversion
    Eh2J = J2Eh ** -1  # Hartree to Joule conversion
    MwVel = []
    if temp is None:
        temp = g09f.getTemp(path_pos)
    AtMass = g09f.getAtMass(path_pos)

    # Extract mass weighted velocities from outputs
    for j in np.nditer(g09f.getVel(path_pos), flags=['external_loop'], order='F'):
        MwVel.append(j)

    Na, Nt, Ns, dt = g09f.getCalcInfo(path_pos)
    t = np.empty(len(MwVel[0]))

    # get time array
    j = 0
    i = 0
    while i < len(t):
        # cartvel[0][i] = np.cos(0.1*i)# for test only
        t[i] = j
        j += dt * 1e-15
        i += 1

    # Smoothing velocities
    # smooth = partial(g09f.smoothFn_windowed, arg=t, Tf=temp, target=2, a=None, b=None)
    # with closing(Pool(processes=8)) as pool:
    #     smvel = pool.map(smooth, MwVel)
    #     pool.terminate()

    # Computing Density of state function
    if calcdos:
        acoriginal = []
        for j in MwVel:
            acoriginal.append(autocorr_manual(j))
        actotone = sumfn(acoriginal)
        dosoriginal = getDos(actotone)

        # ac = []
        # for i in smvel:
        #     ac.append(autocorr_manual(i))
        # acsmvel = sumfn(ac)
        # dossmvel = getDos(acsmvel)
    else:
        pass

    #### Plotting some of the functions if required ####
    if plot and calcdos:
        x1 = range(6000)
        plt.figure(1)
        # plt.subplot(211)
        plt.plot(x1, dosoriginal)
        # plt.subplot(212)
        # plt.plot(x1, dossmvel)
        plt.show()
    else:
        pass

    return dosoriginal  # , dossmvel


def get_Ekinsmavrg_Epotsmavrg(path_pos, path_neg, temp=None, calcHO=False):
    """
    Function to compute average smoothed kinetic and potential energies.
    :param path_pos: Path to positive velocity G09 BOMD out file
    :param path_neg: Path to negative velocity G09 BOMD out file
    :param temp: Temperature of smoothing (default is to get from G09 out files)
    :return: Ep_sm_avr: Smoothed averaged potential energy
             Ek_sm_avr: Smoothed averaged kinetic energy
             T: classical kinetic temperature
    """
    J2Eh = 2.293710449e17  # Joule to Hartree conversion
    Eh2J = J2Eh ** -1  # Hartree to Joule conversion
    MwVelp = []
    MwVeln = []
    MwVel = []
    Forcesp = []
    Forcesn = []
    Forces = []
    coordp = []
    coordn = []
    coord = []
    epot = []
    epotp = []
    epotn = []
    if temp is None:
        temp = g09f.getTemp(path_pos)

    # Extract mass weighted velocities from outputs
    for j in np.nditer(g09f.getVel(path_pos), flags=['external_loop'], order='F'):
        MwVelp.append(j)
    for j in np.nditer(g09f.getVel(path_neg), flags=['external_loop'], order='F'):
        MwVeln.append(j)
    # Extract forces from outputs
    for j in np.nditer(g09f.getF(path_pos), flags=['external_loop'], order='F'):
        Forcesp.append(j)
    for j in np.nditer(g09f.getF(path_neg), flags=['external_loop'], order='F'):
        Forcesn.append(j)
    # Extract potential energy from output
    for j in np.nditer(g09f.getEpot(path_pos), flags=['external_loop'], order='F'):
        epotp.append(j)
    for j in np.nditer(g09f.getEpot(path_neg), flags=['external_loop'], order='F'):
        epotn.append(j)
    # Extract coordinates from output
    for j in np.nditer(g09f.getCoords(path_pos), flags=['external_loop'], order='F'):
        coordp.append(j)
    for j in np.nditer(g09f.getCoords(path_neg), flags=['external_loop'], order='F'):
        coordn.append(j)

    # Concatenate positive and negative arrays
    j = 0
    for i in MwVeln:
        MwVel.append(np.concatenate((np.flipud(i), MwVelp[j][1:])))
        j += 1
    j = 0
    for i in Forcesn:
        Forces.append(np.concatenate((np.flipud(i), Forcesp[j][1:])))
        j += 1
    j = 0
    for i in coordn:
        coord.append(np.concatenate((np.flipud(i), coordp[j][1:])))
        j += 1
    j = 0
    for i in epotn:
        epot.append(np.concatenate((np.flipud(i), epotp[j][1:])))
        j += 1

    Na, Nt, Ns, dt = g09f.getCalcInfo(path_pos)
    t = np.empty(len(MwVel[0]))
    j = 0
    i = 0
    dt = dt * 1e-15
    while i < len(t):
        t[i] = j
        j += dt  # * 1e-15
        i += 1

    smooth = partial(g09f.smoothFn_windowed, arg=t, Tf=temp, target=2, a=None, b=None)
    # Smoothing velocities
    with closing(Pool(processes=8)) as pool:
        smvel = pool.map(smooth, MwVel)
        pool.terminate()
    # Smoothing forces
    with closing(Pool(processes=8)) as pool:
        smforces = pool.map(smooth, Forces)
        pool.terminate()
    # Smoothing coordinates
    with closing(Pool(processes=8)) as pool:
        smcoord = pool.map(smooth, coord)
        pool.terminate()

    # Computing Classical values for heat capacity

    ek = np.transpose(calcEkin(MwVel))
    eksm = np.transpose(calcEkin(smvel))
    # avrEksm = (1 / (len(eksm) * dt)) * integrate.simps(np.multiply(eksm, Eh2J), dx=dt)
    # print(temp,avrEk*2/const.k,avrEksm*2/const.k)
    Ep_sm_avr, Ek_sm_avr, T = getAvrEnergies(epot, eksm, ek, dt, Forces, smforces, coord, smcoord)

    if calcHO:
        smoothHO = partial(g09f.smoothFn_windowed_HO, arg=t, Tf=temp, target=2, a=None, b=None)
        # Smoothing velocities
        with closing(Pool(processes=8)) as pool:
            smvel2 = pool.map(smoothHO, MwVel)
            pool.terminate()
        # Smoothing forces
        with closing(Pool(processes=8)) as pool:
            smforces2 = pool.map(smoothHO, Forces)
            pool.terminate()
        # Smoothing coordinates
        with closing(Pool(processes=8)) as pool:
            smcoord2 = pool.map(smoothHO, coord)
            pool.terminate()

        eksm2 = np.transpose(calcEkin(smvel2))
        Ep_sm_avr2, Ek_sm_avr2, T2 = getAvrEnergies(epot, eksm2, ek, dt, Forces, smforces2, coord, smcoord2)
        return 627503 * Ep_sm_avr, 1000 * Ek_sm_avr, T, 627503 * Ep_sm_avr2, 1000 * Ek_sm_avr2, T2

    # vars = globals().copy()
    # vars.update(locals())
    # shell = code.InteractiveConsole(vars)
    # shell.interact()
    return 627503 * Ep_sm_avr, 1000 * Ek_sm_avr, T


##################################################################
#                           TEST                                 #
##################################################################
def main():
    pass


#    path = 'C:\\Users\\fdavid\\Documents\\entropy\\h2o_bomd\\water_traj99_MD_pos.out'
#    heatcap_velocity(path)

if __name__ == "__main__":
    main()
