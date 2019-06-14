#!/usr/bin/env python3
""""
Here are the utilities for:
 - heat capacity computations
 - Autocorrelation function calculation
 - VDOS calculation
"""

from scipy import linalg, integrate
import scipy.constants as const
import numpy as np
# import modules.g09BOMD_filter as g

R = const.R
NA = const.Avogadro


def autocorr_manual_norm(vel, ts=0.1):
    """
    Function to compute the normalized velocity autocorrelation function
    :param vel: Velocity array
    :param ts: Timstep 
    :return: Normalized velocity autocorrelation function
    """
    ac = np.empty(shape=len(vel), dtype=np.float64)
    ac[0] = 1.0
    norm = integrate.simps(np.multiply(vel, vel), dx=ts * 1E-15)
    for tau in range(1, len(vel)):
        ac[tau] = integrate.simps(np.multiply(vel[tau:], np.roll(vel, tau)[tau:]), dx=ts * 1E-15) / norm
    return ac


def autocorr_manual(vel, ts=0.1):
    """
    Function to compute the velocity autocorrelation function
    :param vel: Velocity array
    :param ts: Timstep 
    :return: Velocity autocorrelation function
    """
    ac = np.empty(shape=len(vel), dtype=np.float64)
    for tau in range(len(vel)):
        ac[tau] = integrate.simps(np.multiply(vel[tau:], np.roll(vel, tau)[tau:]), dx=ts * 1E-15)
    return ac


def sumfn(fn):
    """
    Sums up the given functions for each coordinate
    of each atom
    Args:
        ac: iterable array (list) of numpy arrays containing
         the functions to be summed value by value

   Returns: Numpy array of total autocorrelation function

    """
    tot = np.zeros(shape=len(fn[0]))
    for i in fn:
        k = 0
        for j in i:
            tot[k] = tot[k] + j
            k += 1
    return tot


def getDos(ac, ts=0.1):
    """
    Function to compute Density of state (a.k.a oscillator strength) function from autocorrelation
    Args:
        ac: Autocorrelation function
        ts (float): timestep in fs
    Returns: Density of state as a numpy array in the range of 0-6000 cm-1 with 1 cm-1 step

    """
    dos = np.empty(shape=6000, dtype=np.float64)
    norm = 1 / np.sqrt(2 * const.pi)
    i = 0
    v = 0
    while i < 6000:
        cos = np.cos(-2 * const.pi * v * np.arange(0, len(ac) * 10E-15 * ts, 10E-15 * ts))
        dos[i] = np.multiply(norm, integrate.simps(np.multiply(ac, cos), dx=10E-15 * ts))
        # exp = np.exp(-2 * const.pi * 1j * v * np.arange(0, len(ac) * ts, ts))
        # dos[i] = np.multiply(norm,np.sum(np.multiply(ac,exp)))
        v += 2.99792458E10 * ts
        i += 1
    # dos = np.fft.rfft(ac)
    return dos


def calcEkin(mwvel):
    """
    Calculates kinetic energy from mass weighted! velocities.
    Kinetic energy is simply the square of m.w. velocities
    Args:
        mwvel: iterable array (list) containing numpy arrays of velocity components

    Returns: list of kinetic energy along the trajectory in Hartrees
    """
    # ekin = np.empty(shape=len(mwvel[0]),dtype=np.float64)
    ekin = []
    mwvel = np.transpose(np.array(mwvel))
    for i in mwvel:
        # ek = np.zeros(shape=len(mwvel[0]))
        # for j in i:
        #    ek[k] = j ** 2 * 0.5
        #    k += 1
        # ek = np.divide(np.square(i), 2)
        ek = np.divide(72 * np.square(i), 2 * 9.375828402564380E+29)  # in hartrees
        ekin.append(ek)
        # ek = np.divide(np.sum(np.square(i)),2*9.375828402564380E+29) # in hartrees
        # ekin[l] = ek
        # l += 1

    return ekin


def mw2cart(MwVel, AtMass):
    """
    Transforms mass weighted velocities to cartesian velocities
    Args:
        MwVel: Velocities in mass weighted cartesian coordinates (can be a list of arrays as well)
        AtMass: Array of atomic masses in the same order as they appear
        in MwVel array. If these were extracted by g09BOMD_filter.getVel and
        g09BOMD_filter.getAtMass functions this condition fulfilled by default
    Returns:
        cartvel: velocity in cartesian coordinate system and SI units
    """
    cartvel = []
    for i in MwVel:
        m = 0
        cart = np.empty(shape=len(i))
        k = 0
        for j in i:
            cart[k] = 5.29177249e-11 * np.sqrt(AtMass[m]) * j
            k += 1
            # print(AtMass[m])
            if k % 3 == 0 and k > 0:
                m += 1
                if m == len(AtMass):
                    m = 0
                else:
                    continue
            else:
                continue
        cartvel.append(cart)
    return cartvel


def fitLinear(x, y):
    """
    Function for linear regression 
    :param x: Independent variable array (temperature in this case)
    :param y: Dependent variable array (Total energy)
    :return: Slope, standard deviation 
    """
    if len(x) < 4:
        p = np.polyfit(x, y, deg=1, cov=False)
        return p[0], 0
    p, cov = np.polyfit(x, y, deg=1, cov=True)
    return p[0], np.sqrt(cov[0][0])


def getDos2(ac, ts=0.1):
    """
    Function to compute Density of state (a.k.a oscillator strength) function from autocorrelation
    Args:
        ac: Autocorrelation function
        ts (float): timestep in fs
    Returns: Density of state as a numpy array in the range of 0-6000 cm-1 with 1 cm-1 step
    """
    dos = np.empty(shape=6000, dtype=np.float64)
    norm = 1 / np.sqrt(2 * const.pi)
    i = 0
    v = 0
    while i < 6000:
        cos = np.cos(-2 * const.pi * v * np.arange(0, len(ac) * ts, ts))
        dos[i] = np.multiply(norm, integrate.simps(np.multiply(ac, cos), dx=ts))
        # exp = np.exp(-2 * const.pi * 1j * v * np.arange(0, len(ac) * ts, ts))
        # dos[i] = np.multiply(norm,np.sum(np.multiply(ac,exp)))
        v += 0.1798754748 / (6000.0 * ts)
        i += 1
    # dos = np.fft.rfft(ac)
    return dos


def getDosQM(dos, temp, v0=0.0, ve=6000.0, vstep=1.0):
    """Applies the weighting on the DoS function as proposed by Berens et al."""
    v = np.multiply(2.99792458E10, np.arange(v0, ve, vstep, dtype=np.float64))
    W = np.divide(const.h ** 2 * v ** 2 * np.exp(np.divide(const.h * v, const.k * temp)),
                  const.k ** 2 * temp ** 2 * (1 - np.exp(np.divide(const.h * v, const.k * temp))) ** 2)
    W[0] = 1.0
    dosBerens = np.multiply(dos, W)
    return dosBerens


def getDosG(dos, temp, v0=0.0, ve=6000.0, vstep=1.0):
    """Applies the weighting on the DoS function equivalently to Gaussian smoothing in momentum domain."""
    v = np.multiply(2.99792458E10, np.arange(v0, ve, vstep, dtype=np.float64))
    W = np.exp(np.divide(-const.h ** 2 * v ** 2, 4 * const.pi * const.k ** 2 * temp ** 2))
    dosGauss = np.multiply(dos, W)
    return dosGauss


# def getCV_fromDos(dos, vstep=1.0):
#     """Integrates the DoS over the frequency domain reslting in the heat capacity in J/K*mol"""
#     cv = const.R * integrate.simps(dos, dx=vstep) / 2.99792485E10
#     return cv


def getVirtWork(F, Fsm, x, xsm, dt=1):
    """
    Function to compute virtual work of smoothing. W = 1/2 (F_sm + F) (x_sm-x)
    Args:
        F: Array including force arrays of each atom*3 (x,y,z)
        Fsm: Smoothed force array in the same manner as F
        x: Array including coordinate arrays of each atom*3 (x,y,z)
        xsm: Smoothed coordinates array in the same manner
        dt : Time step of the trajectory points
    Returns: Averaged (1/t * integrate from t0 to t W dt) virtual work for the trajectory
    """

    k = 0
    sum = np.zeros(shape=len(F[0]))
    for fsma in Fsm:
        wi = np.multiply(np.add(fsma, F[k]) * 0.5, (np.add(x[k], np.negative(xsm[k]))))
        sum = np.add(sum, wi)
        k += 1
    w = (1 / ((len(sum) - 1) * dt)) * integrate.simps(sum, dx=dt)
    return w


def getAvrEnergies(Ep, Eksm, Ek, dt, F, Fsm, coord, coordsm):
    """
    Computes average energies and temperature using the virtual work method
    Args:
        Ep: Potential energies array
        Eksm: Smoothed potential energies array
        Ek: Kinetic energies array
        dt: timestep in MD simulations
        F: Forces array
        Fsm: Smoothed forces array
        coord: coordinates array
        coordsm: smoothed coordinates array

    Returns:
        Ep_sm_avr: Smoothed average potential energy computed as:
            <Ep> - W, where <Ep>, and W are the average potential energy 
            and the virtual work of smoothing respectively
        Ek_sm_avr: Smoothed averaged kinetic energy
        T: Classical temperature
    """
    ep_avr = np.float64(0.0)
    eksm_avr = np.float64(0.0)
    ek_avr = np.float64(0.0)
    w = getVirtWork(F, Fsm, coord, coordsm)
    # for ep in Ep[0]:
    #     print(ep)
    for i in Ep:
        ep_avr += (1 / ((len(i) - 1) * dt)) * integrate.simps(i, dx=dt)
    for i in Ek:
        ek_avr += (1 / ((len(i) - 1) * dt)) * integrate.simps(i, dx=dt)
    for i in Eksm:
        eksm_avr += (1 / ((len(i) - 1) * dt)) * integrate.simps(i, dx=dt)
    Ep_sm_avr = np.float64(ep_avr + w)
    Ek_sm_avr = np.float64(eksm_avr * 6.696044921386073e-28)
    T = np.float64(ek_avr * 2 * 2.80162519510793E-24 / ((len(coord)-6) * const.R))
    return Ep_sm_avr, Ek_sm_avr, T

