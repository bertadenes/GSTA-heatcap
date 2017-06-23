#!/usr/bin/env python3
"""
This part controls the calculation of heat capacity.
Using stuff a directory called 'testfiles',  it offers
to change it tough
Starts by reading the kinetic and potential energies from
the given files. Then it calculates the expected kinetic energy
functions on different lower temperatures.
Does the filtering of the kinetic energies, then calculates the
heat capacities both with and without filtering, by differentiation.
"""

from modules.GetHeatcapacity_util import *
from modules.ChangeTemp_util import *
import os
from multiprocessing import Pool
import modules.g09BOMD_filter as g09f
from functools import partial
from contextlib import closing
import warnings
import code
from scipy import integrate
warnings.filterwarnings('ignore')
try:
    import matplotlib.pyplot as plt
    plot = True
except ImportError:
    #print("\nFailed importing matplotlib \nNo plotting, sorry")
    plot = False


# def heatcapacity_processing_2(SubCalc):
#     """Processing classic and smoothed heat capacity calculation
#     Designed for N,V,E (microcanonical ensemble) simulations
#     """
#     # Getting data from VibCalc object :
#     SubCalc.workdir = path
#     #print(path)
#     dT = 10  # Delta temperature for num differentiation
#     time_one = g.time(SubCalc)
#     T = VibCalc.temp + 5
#     Ek_one = g.kineticE(SubCalc)
#     dt = abs(time[1] - time[0])
#     KinE, Ep_period, time = find_period(Ek_one, Ep, time_one)
#     GetAverage(KinE, time)
#     filt_temps = []
#     for i in temps:
#         filt_temps.append(i / 2)
#     j = 0
#     for i in KinEMat:
#         current_time_array = TimeMat[j]
#         print("Start filtering")
#         filt = g.smoothFn(i, current_time_array, filt_temps[j])
#         FiltKinEMat.append(filt)
#         Print2File(filt, "filtered_function" + str(int(temps[j])))
#         j += 1
#
#     AvrTempsfilt, AvrKinEfilt = getAverage(KinE, time)
#
#     Enull = []
#     print("temperatures corresponding to initial kinetic energies are : ", temps)
#     for k in temps:
#         Enull.append(kB * 0.5 * k)
#     print("Average temperatures are : ", AvrTemps)
#     print("The classic heat capacity at ", T, " K is:")
#     CalcCv2(AvrTemps, Enull, T)
#     print("Filt avr  ", AvrTempsfilt)
#     print("The 'filtered' heat capacity at ", T, " K is:")
#     CalcCv2(AvrTempsfilt, Enull, T)


def heatcap_proc_3(path):
    """
    This function calculates the heat capacity using the forces from MD outputs.
    Args:
        path:
    Returns:
    """
    forces = []
    smforces = []
    for j in np.nditer(g09f.GetForces(path), flags=['external_loop'], order ='F'):
        forces.append(j)
    arg = np.arange(len(forces[0]))
    smooth = partial(g09f.smoothFn_ext,arg=arg,Tf=300,target=1,a=None,b=None)
    with closing(Pool(processes=8)) as pool:
        smforces.append(pool.map(smooth, forces))
        pool.terminate()
    for i in smforces:
        for j in i:
            print(j)


def heatcap_velocity(path, calcdos=True, plot=True):
    b2cm = 5.2918e-9
    vel = []
    for j in np.nditer(g09f.getVel(path), flags=['external_loop'], order='F'):
        vel.append(j)
    Na, Nt, Ns, dt = g09f.getCalcInfo(path)
    t = np.empty(len(vel[0]))
    j = 0
    i = 0
    while i < len(t):
        t[i] = j
        j += dt * 1e-15
        i += 1

    # for i in vel:
    #     plt.plot(i)
    #     plt.show()

    if calcdos:
        acoriginal = []
        for j in vel:
            acoriginal.append(autocorr(j))
        actotone = sumfn(acoriginal)
        dosone = getDos(actotone)
    else:
        pass
    smooth = partial(g09f.smoothFn, arg=t, Tf=1500, target=2, a=None, b=None)
    with closing(Pool(processes=8)) as pool:
        # smvel.append(pool.map(smooth, vel))
        smvel = pool.map(smooth, vel)
        pool.terminate()
    if calcdos:
        ac = []
        for i in smvel:
            ac.append(autocorr(i))
        acsum = sumfn(ac)

        # for i in ac:
        #     plt.plot(i)
        #     plt.show()

        dos = getDos(acsum)
    else:
        pass

    #### Plotting some of the functions if required ####
    if plot and calcdos:
        x1 = xaxis(dosone, 33.35641)
        x = xaxis(dos, 33.35641)
        plt.figure(1)
        plt.subplot(211)
        plt.plot(x1, dosone)
        plt.subplot(212)
        plt.plot(x, dos)
        plt.show()
    else:
        pass
    return dos


def heatcap_velocity2(path_pos, path_neg, calcdos=True, plot=True, temp=None):
    J2Eh = 2.293710449e17  # Joule to Hartree conversion
    Eh2J = J2Eh ** -1  # Hartree to Joule conversion
    MwVelp = []
    MwVeln = []
    MwVel = []
    cartvel = []
    Forcesp = []
    Forcesn = []
    Forces = []
    if temp is None:
        temp = g09f.getTemp(path_pos)
    AtMass = g09f.getAtMass(path_pos)

    # Extract mass weighted velocities from outputs
    for j in np.nditer(g09f.getVel(path_pos), flags=['external_loop'], order='F'):
        MwVelp.append(j)
    for j in np.nditer(g09f.getVel(path_neg), flags=['external_loop'], order='F'):
        MwVeln.append(j)

    # Computing cartesian velocities from mass weighted ones
    # cartvelp = mw2cart(MwVelp, AtMass)
    # cartveln = mw2cart(MwVeln, AtMass)

    # SI wrighted velocities
    # cartvelp = []
    # cartveln = []
    # for c in MwVelp:
    #     cartvelp.append(np.sqrt(1.66053904020e-27)*5.29177249e-11*c)
    # for c in MwVeln:
    #     cartveln.append(np.sqrt(1.66053904020e-27)*5.29177249e-11*c)

    #use original velocities, norm later
    cartvelp = MwVelp
    cartveln = MwVeln

    # Extract forces from outputs
    for j in np.nditer(g09f.GetForces(path_pos), flags=['external_loop'], order='F'):
        Forcesp.append(j)
    for j in np.nditer(g09f.GetForces(path_neg), flags=['external_loop'], order='F'):
        Forcesn.append(j)

    # Concatenate positive and negative arrays
    j = 0
    for i in MwVeln:
        MwVel.append(np.concatenate((np.flipud(i), MwVelp[j][1:])))
        j += 1
    j = 0
    for i in cartveln:
        cartvel.append(np.concatenate((np.flipud(i), cartvelp[j][1:])))
        j += 1
    j = 0
    for i in Forcesn:
        Forces.append(np.concatenate((np.flipud(i), Forcesp[j][1:])))
        j +=1

    Na, Nt, Ns, dt = g09f.getCalcInfo(path_pos)
    t = np.empty(len(cartvel[0]))
    # cartvel = []# for test only
    # cartvel.append(np.empty(len(t)))# for test only
    j = 0
    i = 0
    while i < len(t):
        # cartvel[0][i] = np.cos(0.1*i)# for test only
        t[i] = j
        j += dt * 1e-15
        i += 1

    # plt.plot(t, cartvel[0])# for test only
    # plt.show()# for test only

    # Plotting original velocity(time) functions
    # for i in cartvel:
    #     plt.plot(i)
    #     plt.show()

    # Computing Density of state function
    if calcdos:
        acoriginal = []
        for j in cartvel:
            acoriginal.append(autocorr_manual(j))
        actotone = sumfn(acoriginal)
        dosone = getDos(actotone)
    else:
        pass

    # Plotting original velocity autocorrelation functions
    # for i in acoriginal:
    #     plt.plot(i)
    #     plt.show()

    smooth = partial(g09f.smoothFn_windowed, arg=t, Tf=temp, target=2, a=None, b=None)
    # Smoothing velocities
    with closing(Pool(processes=8)) as pool:
        smvel = pool.map(smooth, cartvel)
        pool.terminate()
    # Smoothing forces
    with closing(Pool(processes=8)) as pool:
        smforces = pool.map(smooth, Forces)
        pool.terminate()
    if calcdos:
        ac = []
        for i in smvel:
            ac.append(autocorr_manual(i))
        acsum = sumfn(ac)

        # Plotting smoothed velocity autocorrelation functions
        # for i in ac:
        #     plt.plot(i)
        #     plt.show()

        # plotting the summed autocorrelation functions
        plt.figure(1)
        plt.subplot(211)
        plt.plot(actotone)
        plt.subplot(212)
        plt.plot(acsum)
        plt.show()

        dos = getDos(acsum)
    else:
        pass

    #### Plotting some of the functions if required ####
    if plot and calcdos:
        x1 = range(6000)
        x = xaxis(dosfft, 33.35641)
        plt.figure(1)
        plt.subplot(211)
        plt.plot(x1, dosone)
        plt.subplot(212)
        plt.plot(x1, dos)
        plt.show()
    else:
        pass

    # Computing Classical values for heat capacity
    ek = calcEkin(MwVel)
    ektot = sumfn(ek)

    # Plot some shit
    # plt.plot(ektot)
    # plt.show()


    # EkAvr = getAverage(ektot,t)
    # # amu * bohr^2/ s^2 to Hartree conversion
    # EkAvrEh = np.divide(EkAvr, 9.375828402564380E+29)
    # EkAvrJ = np.multiply(EkAvrEh, Eh2J)

    return dosone, dos

def DOS_from_velocity(path_pos, path_neg, calcdos=True, plot=True, temp=None):
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

    return dosoriginal#, dossmvel

def heatcap_velocity3(SubCalc, calcdos=True, plot=True):
    J2Eh = 2.293710449e17  # Joule to Hartree conversion
    Eh2J = J2Eh ** -1  # Hartree to Joule conversion
    MwVelp = []
    MwVeln = []
    MwVel = []
    cartvel = []
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
    AtMass = g09f.getAtMass(path_pos)

    # Extract mass weighted velocities from outputs
    for j in np.nditer(g09f.getVel(path_pos), flags=['external_loop'], order='F'):
        MwVelp.append(j)
    for j in np.nditer(g09f.getVel(path_neg), flags=['external_loop'], order='F'):
        MwVeln.append(j)

    # Computing cartesian velocities from mass weighted ones
    cartvelp = mw2cart(MwVelp, AtMass)
    cartveln = mw2cart(MwVeln, AtMass)

    # Extract forces from outputs
    for j in np.nditer(g09f.getF(path_pos), flags=['external_loop'], order='F'):
        Forcesp.append(j)
    for j in np.nditer(g09f.getF(path_neg), flags=['external_loop'], order='F'):
        Forcesn.append(j)

    #  Extract coordinates from output
    for j in np.nditer(g09f.getCoords(path_pos), flags=['external_loop'], order='F'):
        coordp.append(j)
    for j in np.nditer(g09f.getCoords(path_neg), flags=['external_loop'], order='F'):
        coordn.append(j)

    #  Extract potential energy from output
    for j in np.nditer(g09f.getEpot(path_pos), flags=['external_loop'], order='F'):
        epotp.append(j)
    for j in np.nditer(g09f.getEpot(path_neg), flags=['external_loop'], order='F'):
        epotn.append(j)

    # Concatenate positive and negative arrays
    j = 0
    for i in MwVeln:
        MwVel.append(np.concatenate((np.flipud(i), MwVelp[j][1:])))
        j += 1
    j = 0
    for i in cartveln:
        cartvel.append(np.concatenate((np.flipud(i), cartvelp[j][1:])))
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
    # From herein the time step is in seconds instead of femtoseconds!
    dt *= 1e-15
    t = np.empty(len(cartvel[0]))
    j = 0
    i = 0
    while i < len(t):
        t[i] = j
        j += dt #* 1e-15
        i += 1

    # Plotting original velocity(time) functions
    # for i in Forces:
    #     plt.plot(i)
    #     plt.show()

    # Computing Density of state function
    if calcdos:
        acoriginal = []
        for j in cartvel:
            acoriginal.append(autocorr_manual(j))
        actotone = sumfn(acoriginal)
        dosone = getDos2(actotone, ts=dt)
    else:
        pass

    # Plotting original velocity autocorrelation functions
    # for i in acoriginal:
    #     plt.plot(i)
    #     plt.show()

    # Smoothing velocities
    smooth = partial(g09f.smoothFn_windowed, arg=t, Tf=temp, target=2, a=None, b=None)
    with closing(Pool(processes=8)) as pool:
        smMwVel = pool.map(smooth, MwVel)
        pool.terminate()

    # Smoothing forces
    with closing(Pool(processes=8)) as pool:
        smforces = pool.map(smooth, Forces)
        pool.terminate()

    # Smoothing coordinates
    # smooth = partial(g09f.smoothFn, arg=t, Tf=temp, target=2, a=None, b=None)
    with closing(Pool(processes=8)) as pool:
        smcoord = pool.map(smooth, coord)
        pool.terminate()

    # Computing kinetic energies

    if calcdos:
        ac = []
        for i in smMwVel:
            ac.append(autocorr_manual(i))
        acsum = sumfn(ac)
        dos = getDos2(acsum, ts=dt)
    else:
        pass

    # Plotting some of the functions if required
    # Plotting smoothed velocity autocorrelation functions
    # for i in ac:
    #     plt.plot(i)
    #     plt.show()
    # plotting the summed autocorrelation functions
    # plt.figure(1)
    # plt.subplot(211)
    # plt.plot(actotone)
    # plt.subplot(212)

    if plot and calcdos:
        plt.figure(1)
        plt.subplot(211)
        plt.plot(dosone)
        plt.subplot(212)
        plt.plot(dos)
        plt.show()
    else:
        pass
    ekin = calcEkin(MwVel)
    smekin = calcEkin(smMwVel)
    return Forces, smforces, coord, smcoord, ekin, smekin, dosone, dos, dt, epot


def GetAvrKinE(ekin, dt):
    ekin_avr = np.float64(0.0)
    for i in ekin:
        ekin_avr += (1/(len(i)*dt))*integrate.simps(i, dx=dt)
    return ekin_avr


def get_Ekinsmavrg_Epotsmavrg(path_pos, path_neg, temp=None):
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
    #  Extract potential energy from output
    for j in np.nditer(g09f.getEpot(path_pos), flags=['external_loop'], order='F'):
        epotp.append(j)
    for j in np.nditer(g09f.getEpot(path_neg), flags=['external_loop'], order='F'):
        epotn.append(j)
    #  Extract coordinates from output
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
        j += dt #* 1e-15
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

    # vars = globals().copy()
    # vars.update(locals())
    # shell = code.InteractiveConsole(vars)
    # shell.interact()
    return 627503*Ep_sm_avr, 1000*Ek_sm_avr, T



##################################################################
#                           TEST                                 #
##################################################################
def main():
    pass
#    path = 'C:\\Users\\fdavid\\Documents\\entropy\\h2o_bomd\\water_traj99_MD_pos.out'
#    heatcap_velocity(path)

if __name__ == "__main__":
    main()