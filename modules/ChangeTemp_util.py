#!/usr/bin/env python3
"""
This program segment converts kinetic energies to different temperatures
Includes some useful functions as well such as Neville interpolation algorithm
and function for exporting data to files.

Update from first version is that it now includes interpolation between certain points
Version 3_0, now returns an equidistant time scale kinetic energy file
Only a module, no program structure

"""
from __future__ import division
from math import *
import scipy.constants as const
import numpy as np

#######################################################
#                   Constants:                        #
#######################################################

kB = const.k  # Boltzmann in J/K
pi = const.pi  # Pi
J2Eh = 2.293710449e17  # Joule to Hartree conversion
Eh2J = J2Eh ** -1  # Hartree to Joule conversion

#######################################################


def T2K(T):
    """Converts kinetic temperature to kinetic energy"""
    return np.float64(1 / 2 * kB * J2Eh * T)


def Neville(vx, vy, x0):
    """
    Efficient algorithm to interpolate an N-1 order polinomial to N points
    Neville(x, y, x0)
    x : array of x values
    y : array of y values
    x0 : you'll get y value for this point
    returns a 64 bit floating point number
    """
    n = len(vx)
    p = n * [0]
    for k in range(n):
        for j in range(n - k):
            if k == 0:
                p[j] = vy[j]
            else:
                p[j] = ((x0 - vx[j + k]) * p[j] + (vx[j] - x0) * p[j + 1]) / (vx[j] - vx[j + k])
    return np.float64(p[0])


def GetTime2(KinE, NewKinE, Time):
    """
    Function to calculate the new time scale from velocities (Kinetic energy arrays)
    GetTime2(KinE, NewKine, Time)
    KinE: kinetic energies on higher temperature
    NewKinE: kinetic energies on lower temperature
    Time: time array corresponding to KinE
    Returns an array
    """
    i = 1
    NewTime = []
    t = 0
    while i < len(NewKinE):
        NewTime.append(t)
        if len(Time) == 0:
            dt = 1
        else:
            dt = Time[i] - Time[i - 1]
        try:
            dtn = np.float64((sqrt(NewKinE[i]) - sqrt(NewKinE[i - 1])) / ((sqrt(KinE[i])) - sqrt(KinE[i - 1])) * dt)
        except ZeroDivisionError:
            dtn = 0.
        t += dtn
        i += 1
    return NewTime

#    i = 1
#    NewTime = np.empty(shape=(len(NewKinE)),dtype=np.float64)#get rid of dynamic allocation
#    t = 0
#    j = 0
#    while i < len(NewKinE):
#        NewTime[j] = t
#        j += 1
#        if len(Time) == 0:
#            dt = 1
#        else:
#            dt = Time[i] - Time[i - 1]
#        try:
#            dtn = np.multiply(np.divide(np.subtract(np.sqrt(NewKinE[i]),np.sqrt(NewKinE[i-1])),np.subtract(np.sqrt(KinE[i]),np.sqrt(KinE[i-1]))),dt)
#        except ZeroDivisionError:
#            dtn = 0.
#        t += dtn
#        i += 1
#    return NewTime


def Print2File(list: object, fname: object) -> object:
    """
    Prints the elements of a list object to a file row by row
    Adds _out to the end of the filename given
    """
    nf = open(fname + "_out.txt", 'w')
    i = 0
    while i < len(list):
        nf.write(str(list[i]) + "\n")
        i = i + 1
    nf.close()
    # print("A new file, named ", fname + "_out", " has been written \n---------------------------")  # at",getcwd())


def Interp2(x, y, y2, dt, EKmax):
    """
    Function to interpolate  between points under 1/4 EkMax, using Neville algorithm
    """
    i = 0
    while i < len(y) - 3:
        step = x[i + 1] - x[i]
        j = 1
        if y[i] <= EKmax / 4:
            while j < int(step / dt):
                p = Neville(x[i:i + 2], y[i:i + 2], x[i] + j * dt)
                q = Neville(x[i:i + 2], y2[i:i + 2], x[i] + j * dt)
                y.insert(i + j, p)
                y2.insert(i + j, q)
                x.insert(i + j, x[i] + j * dt)
                j = j + 1
            i = i + j
        else:
            i = i + 1
    return x, y


def GetEquidistant(Ek, time, dt):
    """Creates new kinetic energy and time arrays now with equidistant time steps"""
    t = 0
    i, j = 0, 0
    NEk = []
    Nt = []
    while i < len(time) - 1:
        if t == time[i]:
            NEk.append(Ek[i])
            Nt.append(time[i])
            t = t + dt
        elif t > time[i] and t < time[i + 1]:
            try:
                ek = Neville(time[i:i + 2], Ek[i:i + 2], t)
                if ek < 0:
                    i = i + 1
                else:
                    NEk.append(ek)
                    Nt.append(t)
                    t = t + dt
            except(ZeroDivisionError):
                i = i + 1
        elif t > time[i + 1]:
            i = i + 1
        else:
            break
    return NEk, Nt


def ReadData():
    """
    Old function to read files needed for the program from argument.
    """
    if len(sys.argv) < 3:
        print(
            "First argument: Kinetic Energy file\nSecond argument: Potential Energy file\nThird argument Time file\nIf no Time file has been given the program assumes equidistant steps")
        exit(1)
    if len(sys.argv) < 4:
        print("Assuming equidistant time steps")
        Time = []
    if len(sys.argv) == 4 or len(sys.argv) == 3:
        try:
            KinE = []
            file = open(sys.argv[1], 'r')
            for line in file:
                KinE.append(float(line))
            file.close()
        except EOFError:
            print("End of file at ", len(data))

        try:
            PotE = []
            file = open(sys.argv[2], 'r')
            for line in file:
                PotE.append(float(line))
            file.close()
        except EOFError:
            print("End of file at ", len(data))
        Time = []
    if len(sys.argv) == 4:
        try:
            Time = []
            file = open(sys.argv[3], 'r')
            for line in file:
                Time.append(float(line))
            file.close()
        except EOFError:
            print("End of file at ", len(data))


def ReadData2(fname):
    """
    Reads data from a  given file row by row and returns with a numpy
    array
    """
    try:
        Data = []
        file = open(fname, 'r')
        for line in file:
            Data.append(float(line))
        file.close()
    except EOFError:
        print("End of file at ", len(data))
    except FileNotFoundError:
        print("Requested file not found!")
        exit(1)
    npData = np.array(Data)
    return npData  # , Data


def MakeTime(arr):
    """
    Creates a time array in the length of the given array filling it
    with whole numbers starting from zero
    """
    j = 0
    time = []
    while j < len(arr):
        time.append(j)
        j += 1
    return time


def MakeTimefs(KinE):
    """
    Creates a time array in the length of the given array filling it
    with numbers starting from zero incrementing by 1e-15
    """
    j = 0
    t = 0
    Time = []
    while j < len(KinE):
        Time.append(t)
        t += 1e-15
        j += 1
    return Time


def CreateTimeArray(dt, len):
    """
    Create an array of given length, and stepsize
    param dt: stepsize
    param len: length of array
    return: numpy array
    """
    i = 0
    time = []
    while i < len:
        time.append(i)
        i += dt
    return np.array(time)


def GetEtot(KinE, PotE):
    """
    Calculates total energy
    param KinE: Kinetic energy array
    param PotE: potential energy array
    return: total energy
    """
    i = 0
    sum = 0.0
    while i < len(KinE):
        y1 = KinE[i]
        y2 = PotE[i]
        sum = sum + y1 + y2
        i = i + 1
    ETot = 1 / (len(KinE)) * sum
    return ETot


def GetKinE(KinE, PotE, T):
    """
    Function to carry out the new kinetic energy calculation on a different temperature
    GetKinE(KinE, PotE, T)
        KinE : Kinetic energy array
        PotE : Potential energy array
        T : New initial(!) temperature *
    Return: New Kinetic energy array, New potential energy array, New value for maximal kinetic energy **
    Use GetTime2() function to calculate the new time array corresponding to the new kinetic energies
         * Calculated as T = 2 * Ek(initial) / kB
         ** Equals to initial kinetic energy
    """
    # print("Please give the new temperature in Kelvin!")
    # T = float(raw_input())
    #print('New average temperature is aimed to be = ', T / 2, " Kelvin")
    K0 = T2K(T)
    #print("Original initial kinetic energy = ", KinE[0], " Hartree")
    #print("New initial kinetic energy = ", K0, " Hartree")
    NewKinE = []
    NewETot = PotE[0] + K0
    #print("New total energy =       ", NewETot, "    Hartree")
    NewPotE = []
    EkMax = max(KinE)
    EkMin = min(KinE)
    NewEkMax = np.float64(K0)
    i = 0
    while i < len(PotE):
        Ki = NewETot - PotE[i]
        if Ki >= 0:
            NewKinE.append(Ki)
            NewPotE.append(PotE[i])
            i += 1
        else:
            NewKinE.append(0.)
            NewPotE.append(PotE[i])
            i += 1
    npnke = np.array(NewKinE)
    npnpe = np.array(NewPotE)
    return NewKinE, NewPotE, NewEkMax


def DriveEkinCalc(KinE, PotE, T, dt, time):
    """
    Manages the Kinetic energy calculation using the functions above
    DriveEkinCalc(KinE, PotE, T, dt, time)
        KinE : Original kinetic energy array (eg. from trajectory)
        PotE : Original potential energy array
        T : Temperature *
        dt = time step for the new kinetic energy array
        time = original time array
    Return: New equidistant kinetic energy array corresponding to the new temperature, New time array
        * Calculated as T = 2 * Ek(initial) / kB
    """
    NewKinE, NewPotE, NewEkMax = GetKinE(KinE, PotE, T)
    NewTime = GetTime2(KinE, NewKinE, time)
    InterpTime, InterpEk = Interp2(NewTime, NewKinE, NewPotE, dt, NewEkMax)
    EqEk, EqT = GetEquidistant(InterpEk, InterpTime, dt)
    return np.array(EqEk), np.array(EqT)


def Go2Zero(traj):
    """
    Modifies a kinetic energy array so it reaches zero at every minimum
    param traj: Trajectory array  (kinetic energies)
    return: numpy array
    """
    i = 1
    while i < len(traj) - 1:
        if traj[i] < traj[i + 1] and traj[i - 1] > traj[i]:
            traj[i] = 0
            i += 1
        else:
            i += 1
    return np.array(traj)

