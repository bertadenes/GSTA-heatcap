#!/usr/bin/env python3
""" version 3.0.0 21/01/2017
    see CHANGE_LOG for details"""
import numpy as np
import scipy.integrate as integrate
import scipy.constants as const
import sys, os, logging, random
from time import sleep
import pickle
from subprocess import Popen
from daemon.g09calcHandler import g09MoveHandler
from daemon.g09calcHandler import g09velgenHandler
from daemon.g09calcHandler import g09MDHandler
import cclib  # parsing gaussian output

from modules.GetHeatcapacity_util import *

try:
    import modules.filelock as fl
except ImportError:
    pass
try:
    import matplotlib.pyplot as plt
except ImportError:
    pass
"""Output process"""


def getCalcInfo(filePath):
    """Read the General parameters section of the output.
    Args:
        filePath (str): The path to the file to process.
    Returns:
        int: Number of atoms in the system.
        int: Number of trajectories in the calculation.
        int: Number of steps taken.
        double: time step if the one used was fixed. Logged if it is so. Otherwise returns 0.
    """
    # version 2

    NAtoms = int(0)
    NTraj = int(0)
    NSteps = int(0)
    dt = np.float64(0)
    with open(filePath, 'r') as outfile:
        for line in outfile:
            if "NAtoms" in line and NAtoms == 0:
                NAtoms = int(line.split()[line.split().index("NAtoms=") + 1])
            elif "Total Number of Trajectories" in line:
                NTraj += int(line.split()[-1])
            elif "Max. points for each Traj." in line and NSteps == 0:
                NSteps += int(line.split()[-1])
            elif "Fixed Trajectory Time Step" in line:
                logging.info("Equidistant time steps used in " + filePath + ":" + line.split()[-2] + "fs.")
                dt = np.float64(line.split()[-2])
    return NAtoms, NTraj, NSteps, dt


def getTrajSummary(filePath):
    """Read Trajectory Summary of g09 BOMD calculation. This section includes the time in each
    step as well as the kinetic and potential energies with smaller accuracy. Returns with an
    array: trajs[trajIndex][lineIndex]=lineString
    Args:
        filePath (str): The path to the file to process.
    Returns:
        list of list of str: Lines read unchanged from the output corresponding to the summary
        for each trajectory.

        First index corresponds to the trajectory, second to a certain line.
        Can be processed by processTrajSummary().
    """
    trajs = []
    i = -1
    read = False
    with open(filePath, 'r') as outfile:
        for line in outfile:
            if "Trajectory summary for trajectory" in line:
                trajs.append([])
                read = True
            elif "******** Start new trajectory calculation ********" in line:
                i += 1
                read = False
            elif "Summary of trajectories" in line:
                read = False
            elif read:
                trajs[i].append(line.strip())
    logging.info("PID" + str(os.getpid()) + ": " + str(len(trajs)) + "trajectories read from file: " + filePath)
    if len(trajs) > 1:
        logging.warn("PID" + str(os.getpid()) + ": More than one trajectory: " + str(len(trajs)))
    return trajs


def processTrajSummary(traj):
    """Processes the summaries for all trajectories read out by getTrajSummary().
    Args:
        traj (list of list of str): Lines read unchanged from the output corresponding to the
        summary for each trajectory.

        First index corresponds to the trajectory, second to a certain line.
        Created by getTrajSummary().
    Returns:
        numpy 2-dimensional array of doubles: Time registerred in the simulation steps in
        femtoseconds.

        First index corresponds to the trajectory, second to the steps.
        numpy 2-dimensional array of doubles: Kinetic energy as given in the summary in a.u.

        First index corresponds to the trajectory, second to the steps.
        WARNING: disables as only given in 1e-7 a.u.
        numpy 2-dimensional array of doubles: Potential energy as given in the summary in a.u.

        First index corresponds to the trajectory, second to the steps.
        WARNING: disables as only given in 1e-7 a.u.
    """
    time = []
    Ekin = []
    Epot = []
    read = False
    for j in range(len(traj)):
        time.append([])
        Ekin.append([])
        Epot.append([])
    for index in range(len(traj)):
        for i in range(len(traj[index])):
            if "Max Error" in traj[index][i]:
                read = False
            if read:
                time[index].append(np.float64(traj[index][i].split()[0]))
                Ekin[index].append(np.float64(traj[index][i].split()[1]))
                Epot[index].append(np.float64(traj[index][i].split()[2]))
            if "Time (fs)   Kinetic (au)   Potent (au)   Delta E (au)     Delta A (h-bar)" in traj[index][
                i]: read = True
    timeArray = np.array(time, dtype=np.float64)
    EkinArray = np.array(Ekin, dtype=np.float64)
    EpotArray = np.array(Epot, dtype=np.float64)
    return timeArray  # ,EkinArray,EpotArray


def getCoords(filePath):
    """Read coordinates in bohr units from g09 BOMD output. Returns with an array like:
    coord[trajIndex][time][atomIndex]=[x,y,z].
    INFO: Redundant as implemented in cclib.
    Args:
        filePath (str): The path to the file to process.
    Returns:
        numpy 4-dimensional array of doubles: Coordinates in angstrom indexed as follows:
        - First index is for trajectories.
        - Second index is for the steps.
        - Third index is for the atom.
        - Fourth index is for the coordinate, 0,1,2 for x,y,z, respectively.
    """
    Na, Nt, Ns, dt = getCalcInfo(filePath)
    coords = np.empty(shape=(Nt, Ns, Na, 3))  # memory allocation
    trajIndex = -1
    readcoord = False
    insum = False
    with open(filePath, 'r') as outfile:
        for line in outfile:
            if "******** Start new trajectory calculation ********" in line:
                trajIndex += 1
                t = -1
            elif "Summary information for step" in line:  # start point and sum 2 are the same
                t += 1
                atomIndex = 0
                insum = True
            elif "Cartesian coordinates: (bohr)" in line:
                readcoord = True
            elif "Predicted information for step" in line or "Normal termination" in line:  # check next step or next link
                insum = False
            elif "I=" not in line:
                readcoord = False
            elif readcoord and insum:
                coords[trajIndex][t][atomIndex][0] = line.split()[3].split('D')[0] + 'E' + line.split()[3].split('D')[1]
                coords[trajIndex][t][atomIndex][1] = line.split()[5].split('D')[0] + 'E' + line.split()[5].split('D')[1]
                coords[trajIndex][t][atomIndex][2] = line.split()[7].split('D')[0] + 'E' + line.split()[7].split('D')[1]
                atomIndex += 1
    # for i in range(len(coords)):
    #     print("Coordinates of",len(coords[i][0]),"atoms were read at",len(coords[i]),"steps for trajectory",i+1,".")

    return coords
    # return 0.52917724900001*coords #get them in angstrom


def getVel(filePath):
    """Read velocities in sqrt(amu)*bohr/sec units from g09 BOMD output. Returns with an array
    like: vel[trajIndex][time][atomIndex]=[v_x,v_y,v_z].
    INFO: Appaerently the atomic velocities are square atomic kinetic energies instead.
    Args:
        filePath (str): The path to the file to process.
    Returns:
        numpy 4-dimensional array of doubles: velocities indexed as follows:
        - First index is for trajectories.
        - Second index is for the steps.
        - Third index is for the atom.
        - Fourth index is for the Cartesian component, 0,1,2 for v_x,v_y,v_z, respectively.
    """
    Na, Nt, Ns, dt = getCalcInfo(filePath)
    vel = np.empty(shape=(Nt, Ns, Na, 3))  # memory allocation
    trajIndex = -1
    readvel = False
    insum = False
    with open(filePath, 'r') as outfile:
        for line in outfile:
            if "******** Start new trajectory calculation ********" in line:
                trajIndex += 1
                t = -1
            elif "Summary information for step" in line:  # start point and sum 2 are the same
                t += 1
                atomIndex = 0
                insum = True
            elif "MW cartesian velocity: (sqrt(amu)*bohr/sec)" in line:
                readvel = True
            elif "Predicted information for step" in line or "Normal termination" in line:  # check next step or next link
                insum = False
            elif "I=" not in line:
                readvel = False
            elif readvel and insum:
                vel[trajIndex][t][atomIndex][0] = line.split()[3].split('D')[0] + 'E' + line.split()[3].split('D')[1]
                vel[trajIndex][t][atomIndex][1] = line.split()[5].split('D')[0] + 'E' + line.split()[5].split('D')[1]
                vel[trajIndex][t][atomIndex][2] = line.split()[7].split('D')[0] + 'E' + line.split()[7].split('D')[1]
                atomIndex += 1
    for i in range(len(vel)):
        logging.info("PID" + str(os.getpid()) + ": Velocities of " + str(len(vel[i][0])) + " atoms were read at " + str(
            len(vel[i])) + " steps for trajectory " + str(i + 1) + ".")
    return vel


def getF(filePath):
    """
    Gets forces at each step of BOMD calculation from G09 output.
    Args:
        filePath (str): The path to the file to process.
    Returns:
        numpy 4-dimensional array of doubles: forces indexed as follows:
        - First index is for trajectories.
        - Second index is for the steps.
        - Third index is for the atom.
        - Fourth index is for the Cartesian component, 0,1,2 for F_x,F_y,F_z, respectively.
    """
    Na, Nt, Ns, dt = getCalcInfo(filePath)
    force = np.empty(shape=(Nt, Ns, Na, 3))  # memory allocation
    trajIndex = 0
    read = False
    t = -1
    finishedTraj1 = False
    with open(filePath, 'r') as outfile:
        for line in outfile:
            if "******** Start new trajectory calculation ********" in line:
                if finishedTraj1:
                    trajIndex += 1
                    t = -1
                else:
                    finishedTraj1 = True
            if "Cartesian Forces:  Max" in line:
                read = False
            if read and len(line.split()) == 5 and line.split()[0] != "Number":
                force[trajIndex][t][atomIndex][0] = np.float64(line.split()[2])
                force[trajIndex][t][atomIndex][1] = np.float64(line.split()[3])
                force[trajIndex][t][atomIndex][2] = np.float64(line.split()[4])
                atomIndex += 1
            if "Forces " in line:
                read = True
                t += 1
                atomIndex = 0
    return force


def getEpot(filePath):
    """Reads the potential energies in hartrees after SCF cycles done.
    Args:
        filePath (str): The path to the file to process.
    Returns:
        numpy 2-dimensional array of doubles: SCF energies indexed as follows:
        - First index is for trajectories.
        - Second index is for the steps.
    """
    Na, Nt, Ns, dt = getCalcInfo(filePath)
    Epot = np.empty(shape=(Nt, Ns), dtype=np.float64)
    trajIndex = -1
    t = -1
    with open(filePath, 'r') as outfile:
        for line in outfile:
            if "******** Start new trajectory calculation ********" in line:
                trajIndex += 1
            elif "SCF Done" in line:
                t += 1
                Epot[trajIndex][t] = line.split()[4]
            elif t == Ns + 1:
                t = -1
    for i in range(len(Epot)):
        logging.info("PID" + str(os.getpid()) + ": SCF potential energies were read at " + str(
            len(Epot[i])) + " steps for trajectory " + str(i + 1) + ".")
    return Epot


def processg09output(filePath, isBOMD=False):
    """Creats a cclib parsed Object with data read from g09 output. Optionally, read BOMD
    related information.
    Args:
        filePath (str): The path to the file to process.
        isBOMD (bool): Should be True if the calculation was BOMD, thus velocities can be read.
    Returns:
        Object: Parsed cclib object. If isBOMD was enabled vel attribute is added to the object
        containing velocity information.
        Also extended with level of theory: data.method and data.basis store the last used level
        of theory.
        Check *http://cclib.github.io/data.html* for further info.
    """
    try:
        try:
            output = cclib.parser.ccopen(filePath)
            data = output.parse()
            with open(filePath, 'r') as f:
                for line in f:
                    if "Standard basis" in line: data.basis = line.split()[2]
                    if "SCF Done" in line: data.method = line.split(')')[0].split('(')[-1]
        except AttributeError:
            print(filePath + " is not a g09 output, or not complete.")
            sys.exit(0)
        if isBOMD:
            data.vel = getVel(filePath)
    except IOError:
        print("File cannot be opened.")
    return data


"""calculation setup"""


def g09setup():
    """Sets up a command for executing Gaussian.
    Returns:
        str: command to submit Gaussian job.
        The command should be recognized by the shell, and its first argument is the input file.
        In addition it is assumed that the output will be directed to the same directory where
        the input was and named the same way as the input was with .log or .out extension.
        It is required for continuing the session after the Gaussian job has terminated.

        Command attribute is stored in the calculation object. At the start of a session, it can
        be given int the config file.
    """
    print("Some information about your Gaussian copy is needed.")
    command = "g09"
    try:
        command = input(
            "Please give the command for executing Gaussian!\nWARNING: output is supposed to be in the same directory "
            "as the input.\n")
    except KeyboardInterrupt:
        print('\n Exiting GSTA')
        sys.exit(0)
    options = ""
    try:
        options = input(
            "Please give the options you wish to specify after the g09 input\nSubmittion will be done: <command> "
            "<input> <options>\n")
    except KeyboardInterrupt:
        print('\n Exiting GSTA')
        sys.exit(0)
    return command, options


def genRndTemp(temp, Natoms):
    """Draw random sample temperature for a Natom-system from chi-square distribution.
    Args:
        temp (np.float64): temperature what must be the expectation value of the distribution
        Natoms (int): the number of atoms in the system.
    Returns:
        rndTemp (np.float64): generated random temperature. The random variable is scaled by T/3N
        in order to have T as expectation value.
    INFO:
        In order to reduce the dependency of the variance of T to the number of atoms, the degrees
        of freedom is given as sqrt(Natoms)


    """
    df = int(np.sqrt(Natoms))
    rndTemp = np.divide(np.random.chisquare(df * 3) * temp, df * 3)
    return rndTemp


def createMoveInputs(Calc):
    """Creats the chosen number of inputs at velgen subdirectory.
    Args:
        Calc (Object): Calculation type object storing the information for a session.
        Only if Calc.type == 2:
            Creates Calc.NTraj number of inputs at Calc.workdir/velgen directory.
    """
    # define elements
    elements = {1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne', 11: 'Na', 12: 'Mg',
                13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca', 21: 'Sc', 22: 'Ti',
                23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu', 30: 'Zn', 31: 'Ga', 32: 'Ge',
                33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr', 37: 'Rb', 38: 'Sr', 39: 'Y', 40: 'Zr', 41: 'Nb', 42: 'Mo',
                43: 'Tc', 44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag', 48: 'Cd', 49: 'In', 50: 'Sn', 51: 'Sb', 52: 'Te',
                53: 'I', 54: 'Xe', 55: 'Cs', 56: 'Ba', 57: 'La', 58: 'Ce', 59: 'Pr', 60: 'Nd', 61: 'Pm', 62: 'Sm',
                63: 'Eu', 64: 'Gd', 65: 'Tb', 66: 'Dy', 67: 'Ho', 68: 'Er', 69: 'Tm', 70: 'Yb', 71: 'Lu', 72: 'Hf',
                73: 'Ta', 74: 'W', 75: 'Re', 76: 'Os', 77: 'Ir', 78: 'Pt', 79: 'Au', 80: 'Hg', 81: 'Tl', 82: 'Pb',
                83: 'Bi', 84: 'Po', 85: 'At', 86: 'Rn', 87: 'Fr', 88: 'Ra', 89: 'Ac', 90: 'Th', 91: 'Pa', 92: 'U',
                93: 'Np', 94: 'Pu', 95: 'Am', 96: 'Cm', 97: 'Bk', 98: 'Cf', 99: 'Es', 100: 'Fm', 101: 'Md', 102: 'No',
                103: 'Lr', 104: 'Rf', 105: 'Db', 106: 'Sg', 107: 'Bh', 108: 'Hs', 109: 'Mt', 110: 'Ds', 111: 'Rg',
                112: 'Cn'}
    dirpath = os.path.join(Calc.workdir, "velgen")
    Calc.setVelGenDir(dirpath)
    if not os.path.exists(dirpath): os.makedirs(dirpath)
    N = Calc.NTraj
    for i in range(N):
        inp = Calc.freqfile.split('.')[0] + "_traj" + str(i + 1) + ".gjf" % (Calc.temp)
        rand = np.random.randint(low=1, high=400000, dtype=np.int64)
        rndtemp = genRndTemp(Calc.temp, Calc.mol.natom)
        com = " BOMD(update=1000,StepSize=1,MaxPoints=1,nsample=" + str(
            len(Calc.mol.vibfreqs)) + ",NTraj=1,Sample=Microcanonical,rtemp=0) IOp(1/44=" + str(rand) + ")"
        title = "Trajectory " + str(i + 1) + " for " + Calc.freqfile.split('.')[0]
        energy = np.float64(const.R * (rndtemp) * 0.000239005736137668)
        Calc.addSubCalc(inp)
        # Calc.SubCalcs[-1].name = "traj"+str(i)
        Calc.SubCalcs[-1].rand = rand
        Calc.SubCalcs[-1].rndtemp = rndtemp
        Calc.SubCalcs[-1].vibFreq = Calc.mol.vibfreqs[0]

        Calc.addMoveInput(inp)

        with open(dirpath + os.sep + inp, 'w') as inputfile:
            inputfile.write("#p " + Calc.mol.method + '/' + Calc.mol.basis + com + "\n\n")
            inputfile.write(title + "\n\n")
            inputfile.write(str(Calc.mol.charge) + ' ' + str(Calc.mol.mult) + '\n')
            for j in range(len(Calc.mol.atomnos)):
                inputfile.write("%s   %f   %f   %f\n" % (
                    elements[Calc.mol.atomnos[j]], Calc.mol.atomcoords[-1][j][0], Calc.mol.atomcoords[-1][j][1],
                    Calc.mol.atomcoords[-1][j][2]))
            inputfile.write("\n1\n")
            for j in range(len(Calc.mol.atomnos)):
                inputfile.write(str(j + 1) + ' ')
            inputfile.write('\n')
            for j in range(len(Calc.mol.vibfreqs)):
                inputfile.write(str(j + 1) + ' ' + str(energy) + ' ')
            inputfile.write('\n\n')
    logging.info("PID" + str(os.getpid()) + ": " + str(N) + " inputs were generated in the " + dirpath + " directory.")
    return


def createMovedVelInputs(SubCalc):
    """
    
    :param SubCalc: 
    :return: 
    """
    # define elements
    elements = {1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne', 11: 'Na', 12: 'Mg',
                13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca', 21: 'Sc', 22: 'Ti',
                23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu', 30: 'Zn', 31: 'Ga', 32: 'Ge',
                33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr', 37: 'Rb', 38: 'Sr', 39: 'Y', 40: 'Zr', 41: 'Nb', 42: 'Mo',
                43: 'Tc', 44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag', 48: 'Cd', 49: 'In', 50: 'Sn', 51: 'Sb', 52: 'Te',
                53: 'I', 54: 'Xe', 55: 'Cs', 56: 'Ba', 57: 'La', 58: 'Ce', 59: 'Pr', 60: 'Nd', 61: 'Pm', 62: 'Sm',
                63: 'Eu', 64: 'Gd', 65: 'Tb', 66: 'Dy', 67: 'Ho', 68: 'Er', 69: 'Tm', 70: 'Yb', 71: 'Lu', 72: 'Hf',
                73: 'Ta', 74: 'W', 75: 'Re', 76: 'Os', 77: 'Ir', 78: 'Pt', 79: 'Au', 80: 'Hg', 81: 'Tl', 82: 'Pb',
                83: 'Bi', 84: 'Po', 85: 'At', 86: 'Rn', 87: 'Fr', 88: 'Ra', 89: 'Ac', 90: 'Th', 91: 'Pa', 92: 'U',
                93: 'Np', 94: 'Pu', 95: 'Am', 96: 'Cm', 97: 'Bk', 98: 'Cf', 99: 'Es', 100: 'Fm', 101: 'Md', 102: 'No',
                103: 'Lr', 104: 'Rf', 105: 'Db', 106: 'Sg', 107: 'Bh', 108: 'Hs', 109: 'Mt', 110: 'Ds', 111: 'Rg',
                112: 'Cn'}
    dirpath = SubCalc.workdir
    inp = SubCalc.freqfile.split('.')[0] + '_' + SubCalc.name + "_velgen.gjf"
    SubCalc.velgenInput = inp
    rtemp = '0'
    SubCalc.rtemp = int(rtemp)
    com = " BOMD(update=1000,StepSize=1,MaxPoints=2,Sample=Fixed,NSample=" + str(
        len(SubCalc.mol.modeEkin)) + ",NTraj=1,rtemp=" + rtemp + ") nosymmetry"
    title = SubCalc.name + " for " + SubCalc.freqfile.split('.')[0]
    with open(dirpath + os.sep + inp, 'w') as inputfile:
        inputfile.write("#p " + SubCalc.mol.method + '/' + SubCalc.mol.basis + com + "\n\n")
        inputfile.write(title + "\n\n")
        inputfile.write(str(SubCalc.mol.charge) + ' ' + str(SubCalc.mol.mult) + '\n')
        for j in range(len(SubCalc.mol.atomnos)):
            inputfile.write("%s   %f   %f   %f\n" % (
                elements[SubCalc.mol.atomnos[j]], SubCalc.mol.atomcoords[-1][j][0], SubCalc.mol.atomcoords[-1][j][1],
                SubCalc.mol.atomcoords[-1][j][2]))
        inputfile.write("\n1\n")
        for j in range(len(SubCalc.mol.atomnos)):
            inputfile.write(str(j + 1) + ' ')
        inputfile.write('\n')
        for j in range(len(SubCalc.mol.modeEkin)):
            inputfile.write(str(j + 1) + ' ' + str(SubCalc.mol.modeEkin[j]) + ' ')
        inputfile.write('\n\n')
    logging.info(
        "PID" + str(os.getpid()) + ": one velgen input was generated in the " + dirpath + " directory for " + str(
            SubCalc.name))
    return


def createVelInputs(Calc):
    """Creates the chosen number of inputs at velgen subdirectory.
    Args:
        Calc (Object): Calculation type object storing the information for a session.

        If Calc.type == 1:
            Creates Calc.maxvib number of inputs at Calc.workdir/velgen directory. For each
            calculation, a VibCalc object is constructed.

        If Calc.type == 2:
            Creates Calc.NTraj number of inputs at Calc.workdir/velgen directory.
    INFO: for vibwise calculation the energy on the modes chosen according to the temperature
    given, however it is 5 K higher. It is due to the fact that PES at x K can be derived from
    PES at x+5 K, but not the other way around.
    It cannot be done for NVE sampling thus 3 inputs created (T-5,T,T+5) with same random seed.
    """
    # define elements
    elements = {1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne', 11: 'Na', 12: 'Mg',
                13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca', 21: 'Sc', 22: 'Ti',
                23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu', 30: 'Zn', 31: 'Ga', 32: 'Ge',
                33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr', 37: 'Rb', 38: 'Sr', 39: 'Y', 40: 'Zr', 41: 'Nb', 42: 'Mo',
                43: 'Tc', 44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag', 48: 'Cd', 49: 'In', 50: 'Sn', 51: 'Sb', 52: 'Te',
                53: 'I', 54: 'Xe', 55: 'Cs', 56: 'Ba', 57: 'La', 58: 'Ce', 59: 'Pr', 60: 'Nd', 61: 'Pm', 62: 'Sm',
                63: 'Eu', 64: 'Gd', 65: 'Tb', 66: 'Dy', 67: 'Ho', 68: 'Er', 69: 'Tm', 70: 'Yb', 71: 'Lu', 72: 'Hf',
                73: 'Ta', 74: 'W', 75: 'Re', 76: 'Os', 77: 'Ir', 78: 'Pt', 79: 'Au', 80: 'Hg', 81: 'Tl', 82: 'Pb',
                83: 'Bi', 84: 'Po', 85: 'At', 86: 'Rn', 87: 'Fr', 88: 'Ra', 89: 'Ac', 90: 'Th', 91: 'Pa', 92: 'U',
                93: 'Np', 94: 'Pu', 95: 'Am', 96: 'Cm', 97: 'Bk', 98: 'Cf', 99: 'Es', 100: 'Fm', 101: 'Md', 102: 'No',
                103: 'Lr', 104: 'Rf', 105: 'Db', 106: 'Sg', 107: 'Bh', 108: 'Hs', 109: 'Mt', 110: 'Ds', 111: 'Rg',
                112: 'Cn'}
    dirpath = os.path.join(Calc.workdir, "velgen")
    Calc.setVelGenDir(dirpath)
    if not os.path.exists(dirpath): os.makedirs(dirpath)
    N = Calc.NTraj
    for i in range(N):
        inp = Calc.freqfile.split('.')[0] + "_traj" + str(i + 1) + ".gjf" % (Calc.temp)
        rand = np.random.randint(low=1, high=400000, dtype=np.int64)
        rndtemp = genRndTemp(Calc.temp, Calc.mol.natom)
        rtemp = str(int(rndtemp))
        com = " BOMD(update=1000,StepSize=1,MaxPoints=1,nsample=" + str(
            len(Calc.mol.vibfreqs)) + ",NTraj=1,Sample=Microcanonical,rtemp=" + rtemp + ") IOp(1/44=" + str(rand) + ")"
        title = "Trajectory " + str(i + 1) + " for " + Calc.freqfile.split('.')[0]
        energy = np.float64(int(rndtemp) * 0.001987203611)
        Calc.addSubCalc(inp)
        Calc.SubCalcs[-1].rand = rand
        Calc.SubCalcs[-1].rtemp = int(rndtemp)
        Calc.SubCalcs[-1].rndtemp = rndtemp
        Calc.SubCalcs[-1].vibFreq = Calc.mol.vibfreqs[0]

        Calc.addVelGenInput(inp)

        with open(dirpath + os.sep + inp, 'w') as inputfile:
            inputfile.write("%chk=" + inp.split('.')[0] + ".chk\n\n")
            inputfile.write("#p " + Calc.mol.method + '/' + Calc.mol.basis + com + "\n\n")
            inputfile.write(title + "\n\n")
            inputfile.write(str(Calc.mol.charge) + ' ' + str(Calc.mol.mult) + '\n')
            for j in range(len(Calc.mol.atomnos)):
                inputfile.write("%s   %f   %f   %f\n" % (
                    elements[Calc.mol.atomnos[j]], Calc.mol.atomcoords[-1][j][0], Calc.mol.atomcoords[-1][j][1],
                    Calc.mol.atomcoords[-1][j][2]))
            inputfile.write("\n1\n")
            for j in range(len(Calc.mol.atomnos)):
                inputfile.write(str(j + 1) + ' ')
            inputfile.write('\n')
            for j in range(len(Calc.mol.vibfreqs)):
                inputfile.write(str(j + 1) + ' ' + str(energy) + ' ')
            inputfile.write("\n\n")

    logging.info("PID" + str(os.getpid()) + ": " + str(N) + " inputs were generated in the " + dirpath + " directory.")
    return


def runVelGen(Calc):
    if os.fork():
        logging.info("PID" + str(os.getpid()) + ": Exiting entropy\n")
        sys.exit(0)
    VH = g09velgenHandler(Calc)
    VH.watch()
    return


def runMove(Calc):
    if os.fork():
        logging.info("PID" + str(os.getpid()) + ": Exiting entropy\n")
        sys.exit(0)
    MH = g09MoveHandler(Calc)
    MH.watch()
    return


def createMDInput(subCalc):
    """Creates pair of inputs for BOMD simulation either along one certain vibration in both
    direction or for one independent trajectory.
    SubCalc object must resemble velocities, it can be updated using updateSubCalc().
    Args:
        subCalc (Object): SubCalc type object storing information for one vibration or trajectory.

        Creates 2 inputs in subCalc.workdir ending _pos and _neg, corresponding if the velocities
        were multiplied by -1 or not. subCalc.mol must have an attribute vel, which is added by
        calling processg09output() with isBOMD=True option. The names of the inputs are stored
        in subCalc.MDinput1 and subCalc.MDinput2.
    """
    elements = {1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne', 11: 'Na', 12: 'Mg',
                13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca', 21: 'Sc', 22: 'Ti',
                23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu', 30: 'Zn', 31: 'Ga', 32: 'Ge',
                33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr', 37: 'Rb', 38: 'Sr', 39: 'Y', 40: 'Zr', 41: 'Nb', 42: 'Mo',
                43: 'Tc', 44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag', 48: 'Cd', 49: 'In', 50: 'Sn', 51: 'Sb', 52: 'Te',
                53: 'I', 54: 'Xe', 55: 'Cs', 56: 'Ba', 57: 'La', 58: 'Ce', 59: 'Pr', 60: 'Nd', 61: 'Pm', 62: 'Sm',
                63: 'Eu', 64: 'Gd', 65: 'Tb', 66: 'Dy', 67: 'Ho', 68: 'Er', 69: 'Tm', 70: 'Yb', 71: 'Lu', 72: 'Hf',
                73: 'Ta', 74: 'W', 75: 'Re', 76: 'Os', 77: 'Ir', 78: 'Pt', 79: 'Au', 80: 'Hg', 81: 'Tl', 82: 'Pb',
                83: 'Bi', 84: 'Po', 85: 'At', 86: 'Rn', 87: 'Fr', 88: 'Ra', 89: 'Ac', 90: 'Th', 91: 'Pa', 92: 'U',
                93: 'Np', 94: 'Pu', 95: 'Am', 96: 'Cm', 97: 'Bk', 98: 'Cf', 99: 'Es', 100: 'Fm', 101: 'Md', 102: 'No',
                103: 'Lr', 104: 'Rf', 105: 'Db', 106: 'Sg', 107: 'Bh', 108: 'Hs', 109: 'Mt', 110: 'Ds', 111: 'Rg',
                112: 'Cn'}
    title = "NVE sampling seed " + str(subCalc.rand) + " T = " + str(int(subCalc.rndtemp))
    stepsize = int(10E14 / (subCalc.vibFreq * 1199169832))

    # set to 0.1 fs
    stepsize = 1000
    d = getDecay(subCalc.temp)
    MP = str(int(d / (stepsize * 10E-20)) + 201)
    # set to sum up to 1 ps
    MP = str(5000)
    subCalc.MP = MP
    if subCalc.rot:
        rtemp = str(int(subCalc.rndtemp))
    else:
        rtemp = str(0)
    dirpath = subCalc.workdir
    if not os.path.exists(dirpath): os.makedirs(dirpath)
    with open(dirpath + os.sep + subCalc.inputname.split('.')[0] + "_MD_pos.gjf", 'w') as inputfile:
        subCalc.MDinput1 = subCalc.inputname.split('.')[0] + "_MD_pos.gjf"
        inputfile.write("#p " + subCalc.mol.method + '/' + subCalc.mol.basis + " BOMD(gradientonly,StepSize=" + str(
            stepsize) + ",MaxPoints=" + MP + ",ReadMWVelocity,rtemp=" + rtemp + ") nosymmetry\n\n")
        inputfile.write(title + "\n\n")
        inputfile.write(str(subCalc.mol.charge) + ' ' + str(subCalc.mol.mult) + '\n')
        for j in range(len(subCalc.mol.atomnos)):
            inputfile.write("%s   %f   %f   %f\n" % (
                elements[subCalc.mol.atomnos[j]], subCalc.mol.atomcoords[-1][j][0], subCalc.mol.atomcoords[-1][j][1],
                subCalc.mol.atomcoords[-1][j][2]))
        inputfile.write("\n1\n")
        for j in range(len(subCalc.mol.atomnos)):
            inputfile.write(str(j + 1) + ' ')
        inputfile.write('\n')
        for j in range(len(subCalc.mol.atomnos)):
            inputfile.write(
                str(subCalc.mol.vel[-1][-1][j][0]) + "   " + str(subCalc.mol.vel[-1][-1][j][1]) + "   " + str(
                    subCalc.mol.vel[-1][-1][j][2]) + "\n")
        inputfile.write('\n')
    with open(dirpath + os.sep + subCalc.inputname.split('.')[0] + "_MD_neg.gjf", 'w') as inputfile:
        subCalc.MDinput2 = subCalc.inputname.split('.')[0] + "_MD_neg.gjf"
        inputfile.write("#p " + subCalc.mol.method + '/' + subCalc.mol.basis + " BOMD(gradientonly,StepSize=" + str(
            stepsize) + ",MaxPoints=" + MP + ",ReadMWVelocity,rtemp=" + rtemp + ") nosymmetry\n\n")
        inputfile.write(title + "\n\n")
        inputfile.write(str(subCalc.mol.charge) + ' ' + str(subCalc.mol.mult) + '\n')
        for j in range(len(subCalc.mol.atomnos)):
            inputfile.write("%s   %f   %f   %f\n" % (
                elements[subCalc.mol.atomnos[j]], subCalc.mol.atomcoords[-1][j][0], subCalc.mol.atomcoords[-1][j][1],
                subCalc.mol.atomcoords[-1][j][2]))
        inputfile.write("\n1\n")
        for j in range(len(subCalc.mol.atomnos)):
            inputfile.write(str(j + 1) + ' ')
        inputfile.write('\n')
        for j in range(len(subCalc.mol.atomnos)):
            inputfile.write(
                str(-1 * subCalc.mol.vel[-1][-1][j][0]) + "   " + str(-1 * subCalc.mol.vel[-1][-1][j][1]) + "   " + str(
                    -1 * subCalc.mol.vel[-1][-1][j][2]) + "\n")
        inputfile.write('\n')
    logging.info("PID" + str(os.getpid()) + ": 2 inputs were generated in " + dirpath + " for " + subCalc.name + ".")
    return


def runMD(Calc: object) -> object:
    if os.fork():
        logging.info("PID" + str(os.getpid()) + ": Exiting entropy\n")
        sys.exit(0)
    MH = g09MDHandler(Calc)
    MH.watch()
    return


"""data process"""


def time(SC):
    """Gets the time from a pair of BOMD calculations. As the two calculations had identical
    number of steps and stepsize, but in opposite direction, the time range will be (-t,t) with
    equidistant steps. Thus if the simulation was N step long, the array will be 2N-1 long
    because 0 presents twice.
    Args:
        SC (Object): SubCalc type object storing the information for one vibration or trajectory.

        SC must have MDoutput1 and MDoutput2 attributes, which can be added calling
        updateMDpos() and updateMDneg() functions respectively.
    Returns:
        numpy array of doubles: time values in femtoseconds in a range of (-t,t)

    WARNING: it only works with calculations resembling 1 trajectory.
    """
    traj_pos = getTrajSummary(SC.MDoutput1)
    traj_neg = getTrajSummary(SC.MDoutput2)
    t_pos = processTrajSummary(traj_pos)
    t_neg = processTrajSummary(traj_neg)
    t = np.concatenate((-t_neg[0][:0:-1], t_pos[0]))  # WARNING: only working with 1 trajectory
    return t


def kineticE(SC):
    """Calculates and returns the kinetic energy using corrected atomic velocities. Components
    are projected to the initial velocity vector filtering dissipated energy.
    Args:
        SC (Object): VibCalc type object storing the information for one vibration or trajectory.

        SC must have MDoutput1 and MDoutput2 attributes, which can be added calling
        updateMDpos() and updateMDneg() functions respectively. Then the data of the second
        output is reversed in time counting down from 0.
    Returns:
        numpy array of doubles: kinetic energies in hartrees respect to one normal mode,
        derived from the atomic velocities. 

    WARNING: it only works with calculations consisting 1 trajectory.
    """
    os.chdir(SC.workdir)
    v_pos = getVel(SC.MDoutput1)
    v_neg = getVel(SC.MDoutput2)
    v0_len = np.linalg.norm(v_pos[0][0][:][:], axis=1)  # WARNING: only working with 1 trajectory
    v0_norm = np.divide(v_pos[0][0][:][:], v0_len.reshape(-1, 1))
    E_kinpos_atom = np.empty(shape=(len(v_pos[0]), len(v_pos[0][0])), dtype=np.float64)
    E_kinneg_atom = np.empty(shape=(len(v_neg[0]), len(v_neg[0][0])), dtype=np.float64)
    i = 0
    for t in v_pos[0][:]:
        j = 0
        for a in t:
            E_kinpos_atom[i][j] = np.divide(np.power(np.dot(v_pos[0][i][j], v0_norm[j]), 2), 2)
            j += 1
        i += 1
    i = 0
    for t in v_neg[0][:]:
        j = 0
        for a in t:
            E_kinneg_atom[i][j] = np.divide(np.power(np.dot(v_neg[0][i][j], v0_norm[j]), 2), 2)
            j += 1
        i += 1
    E_kinpos = np.sum(E_kinpos_atom, axis=1)
    E_kinneg = np.sum(E_kinneg_atom, axis=1)
    E_kin = np.divide(np.concatenate((E_kinneg[:0:-1], E_kinpos)), 9.375828402564380E+29)
    return E_kin


def getDecay(Tf, b=None):
    kB = const.k
    h = const.h
    pi = const.pi
    B = np.float64((-8 * kB ** 2 * Tf ** 2 * pi ** 3) / (h ** 2))
    # get gaussian decay
    d = np.sqrt(np.divide(np.log(np.float64(1e-8)), B))
    return d


def getAtMass(path):
    """
    This function extracts the atomic mass data from a Gaussian output file.
    Returns an array of atomic masses for each atom in the same order as it
    appeared in the output.
    Args:
        path: Path to the output file

    Returns:
        AtMass: Array of atomic masses

    """
    file = open(path, "r")
    Na, Nt, Ns, dt = getCalcInfo(path)
    AtMass = np.empty(shape=Na)
    i = 0
    for line in file:
        if "AtmWgt=" in line:
            for k in line.split()[1:]:
                AtMass[i] = k
                if i == Na - 1:
                    break
                else:
                    i += 1
        else:
            continue
    return AtMass


def getTemp(path):
    file = open(path, "r")
    Na, Nt, Ns, dt = getCalcInfo(path)
    AtMass = np.empty(shape=Na)
    T = np.float64(0)
    i = 0
    for line in file:
        if "NVE sampling seed" in line:
            T = line.split()[6]
            break
    return float(T)


class NonSmoothableTraj(Exception):
    pass


def smoothFn_fixstep(fn, step, Tf, target=1, a=None, b=None):
    kB = const.k
    h = const.h
    pi = const.pi
    if target == 1:
        # for energies
        A = np.float64((np.sqrt(32) * pi * kB * Tf) / h)
        B = np.float64((-32 * kB ** 2 * Tf ** 2 * pi ** 3) / (h ** 2))
    elif target == 2:
        # for coordinates and velocities
        A = np.float64((np.sqrt(8) * pi * kB * Tf) / h)
        B = np.float64((-8 * kB ** 2 * Tf ** 2 * pi ** 3) / (h ** 2))
    else:
        A, B = a, b
    smfn = arg = np.empty(shape=len(fn), dtype=np.float64)
    for i in range(len(arg)):
        arg[i] = i * step
    G = np.multiply(A, np.exp(np.multiply(np.square(arg), B)))
    smfn[0] = integrate.simps(np.multiply(fn, G), dx=step) / integrate.simps(G, dx=step)
    for i in range(1, len(fn)):
        smfn[i] = integrate.simps(np.multiply(fn[i:], G[:-i]), dx=step) / integrate.simps(G[:-i],
                                                                                          dx=step) + integrate.simps(
            np.multiply(fn[:i], np.flipud(G)[-i - 1:-1]), dx=step) / integrate.simps(np.flipud(G)[-i - 1:-1], dx=step)
    return smfn


def smoothFn(fn, arg, Tf, target=1, a=None, b=None):
    """Filter a function along one argument using a Gaussian function *g* with adjustable
    parameters *A* and *B*.
    `g(tau) = A * exp(B * tau^2)`

    It uses the whole set of data points for smoothing but only returns a truncated set, where
    of Gaussian decay left from both sides.

    Args:
        fn (numpy array of doubles): the values of the function at discrete points.
        arg (numpy array of doubles): the argument values corresponding to the function values.
        Smoothing will be done with respect to this argument.
        Tf (double): temperature for filtering.
        target (int): specifies the nature of the function to be smoothed.
        1 stands for energies, 2 for coordinates or velocities, needs to be manually given
        otherwise.
        a (double): gives the value of *A* if manually given.
        b (double): gives the value of *B* if manually given.
    Returns:
        numpy array of doubles: smoothed function at the same points where the original was
        given.
    """
    kB = const.k
    h = const.h
    pi = const.pi
    if target == 1:
        # for energies
        A = np.float64((np.sqrt(32) * pi * kB * Tf) / h)
        B = np.float64((-32 * kB ** 2 * Tf ** 2 * pi ** 3) / (h ** 2))
    elif target == 2:
        # for coordinates and velocities
        A = np.float64((np.sqrt(8) * pi * kB * Tf) / h)
        B = np.float64((-8 * kB ** 2 * Tf ** 2 * pi ** 3) / (h ** 2))
    else:
        A, B = a, b

    # get gaussian decay
    d = np.sqrt(np.divide(np.log(np.float64(1e-8)), B))
    ext = int(d / (arg[1] - arg[0])) + 1

    # smarg = arg[ext:-ext]
    # smarg = np.empty(shape=len(arg)-2*ext,dtype=np.float64)
    try:
        smfn = np.empty(shape=len(arg) - 2 * ext, dtype=np.float64)
    except ValueError:
        raise NonSmoothableTraj
    for i in range(len(arg)):
        G = np.multiply(A, np.exp(np.multiply(np.square(np.subtract(arg, arg[i])), B)))
        fng = np.multiply(fn, G)
        if i > ext - 1 and i < len(arg) - ext:
            # smarg[i-ext] = arg[i]
            smfn[i - ext] = integrate.simps(fng, dx=arg[1] - arg[0])
    # return smarg,smfn
    return smfn


def smoothFn_ext(fn, arg, Tf, target=1, a=None, b=None):
    """Filter a function along one argument using a Gaussian function *g* with adjustable
    parameters *A* and *B*.
    `g(tau) = A * exp(B * tau^2)`

    During smoothing process, the function is extended by a period in which the smoothing
    function decays significally (below 1E-10*maximum). This is to avoid artifacts in the
    beginning and the end of the trajectory.

    Args:
        fn (numpy array of doubles): the values of the function at discrete points.
        arg (numpy array of doubles): the argument values corresponding to the function values.
        Smoothing will be done with respect to this argument.
        Tf (double): temperature for filtering.
        target (int): specifies the nature of the function to be smoothed.
        1 stands for energies, 2 for coordinates or velocities, needs to be manually given
        otherwise.
        a (double): gives the value of *A* if manually given.
        b (double): gives the value of *B* if manually given.
    Returns:
        numpy array of doubles: smoothed function at the same points where the original was
        given.

    """
    kB = const.k
    h = const.h
    pi = const.pi
    if target == 1:
        # for energies
        A = np.float64((np.sqrt(32) * pi * kB * Tf) / h)
        B = np.float64((-32 * kB ** 2 * Tf ** 2 * pi ** 3) / (h ** 2))
    elif target == 2:
        # for coordinates and velocities
        A = np.float64((np.sqrt(8) * pi * kB * Tf) / h)
        B = np.float64((-8 * kB ** 2 * Tf ** 2 * pi ** 3) / (h ** 2))
    else:
        A, B = a, b

    # get gaussian decay
    d = np.divide(np.sqrt(np.divide(np.log(np.float64(1e-10)), B)), 2)
    ext = int(d / (arg[1] - arg[0])) + 1
    # extending arg by ext in both direction
    argEx = np.empty(shape=len(arg) + 2 * ext, dtype=np.float64)
    for i in range(ext):
        argEx[i] = arg[0] - (arg[-1] - arg[(-ext + i - 1) % len(arg)])
    for i in range(ext, len(arg) + ext):
        argEx[i] = arg[i - ext]
    for i in range(len(arg) + ext, len(arg) + 2 * ext):
        argEx[i] = arg[-1] + (arg[(i + 1 - len(arg) - ext) % len(arg)] - arg[0])
    # extending fn  by ext in both direction
    fnEx = np.empty(shape=len(arg) + 2 * ext, dtype=np.float64)
    for i in range(len(fnEx)):
        fnEx[i] = fn[i % len(arg) - ext + 1]

    smfn = np.empty(shape=len(arg), dtype=np.float64)
    k = 0
    for i in range(len(argEx)):
        G = np.multiply(A, np.exp(np.multiply(np.square(np.subtract(argEx, argEx[i])), B)))
        fng = np.multiply(fnEx, G)
        integral = np.float64(0.0)
        if i > ext - 1 and i < len(arg) + ext:
            for j in range(len(argEx) - 1):
                # (fng[j+1]+fng[j])*(argEx[j+1]-argEx[j])/2
                integral += np.divide(np.multiply(np.add(fng[j + 1], fng[j]), np.subtract(argEx[j + 1], argEx[j])), 2)
            smfn[k] = integral
            k += 1
    return smfn


def smoothFn_old(fn, arg, Tf, target=1, a=None, b=None):
    """g(tau) = A * exp(B * tau^2)"""
    kB = const.k
    h = const.h
    pi = const.pi
    # Tf = np.float64(300)
    # try:
    #    Tf = np.float64(input("Please give the temperature for the filtering (Tf)!\n"))
    # except ValueError:
    #    print("Please give a valid number!")
    # Tf read rather as arg for multimple function calls
    if target == 1:
        # for energies
        A = np.float64((np.sqrt(32) * pi * kB * Tf) / h)
        B = np.float64((-32 * kB ** 2 * Tf ** 2 * pi ** 3) / (h ** 2))
    elif target == 2:
        # for coordinates and velocities
        A = np.float64((np.sqrt(8) * pi * kB * Tf) / h)
        B = np.float64((-8 * kB ** 2 * Tf ** 2 * pi ** 3) / (h ** 2))
    else:
        A, B = a, b
    smfn = np.empty(shape=len(arg), dtype=np.float64)
    for i in range(len(arg)):
        G = np.multiply(A, np.exp(np.multiply(np.square(np.subtract(arg, arg[i])), B)))
        norm = np.sum(G) * (arg[1] - arg[0])
        fng = np.divide(np.multiply(fn, G), norm)
        integral = np.float64(0.0)
        for j in range(len(arg) - 1):
            # (fng[j+1]+fng[j])*(arg[j+1]-arg[j])/2
            integral += np.divide(np.multiply(np.add(fng[j + 1], fng[j]), np.subtract(arg[j + 1], arg[j])), 2)
        smfn[i] = integral
    return smfn


def smoothFn_windowed(fn, arg, Tf, target=2, a=None, b=None):
    """
    formula: g(tau) = A * exp(B * tau^2)
    :param fn: 
    :param arg: 
    :param Tf: 
    :param target: 
    :param a: 
    :param b: 
    :return: 
    """
    kB = const.k
    h = const.h
    pi = const.pi
    # Tf = np.float64(300)
    # try:
    #    Tf = np.float64(input("Please give the temperature for the filtering (Tf)!\n"))
    # except ValueError:
    #    print("Please give a valid number!")
    # Tf read rather as arg for multimple function calls
    if target == 1:
        # for energies
        A = np.float64((np.sqrt(32) * pi * kB * Tf) / h)
        B = np.float64((-32 * kB ** 2 * Tf ** 2 * pi ** 3) / (h ** 2))
    elif target == 2:
        # for coordinates and velocities
        A = np.float64((np.sqrt(8) * pi * kB * Tf) / h)
        B = np.float64((-8 * kB ** 2 * Tf ** 2 * pi ** 3) / (h ** 2))
    else:
        A, B = a, b
    # numfile = open("numfn10e8.bin", "rb")
    # sm = pickle.load(numfile)
    # numfile.close()

    smfn = np.empty(shape=len(arg), dtype=np.float64)
    for i in range(len(arg)):
        G = np.multiply(A, np.exp(np.multiply(np.square(np.subtract(arg, arg[i])), B)))
        # if i == 0:
        #     print(G)
        #     plt.plot(arg, sm, 'x', arg, G, "o")
        #     plt.show()
        norm = integrate.simps(G, dx=arg[1] - arg[0])
        fng = np.multiply(fn, G)
        smfn[i] = integrate.simps(fng, dx=arg[1] - arg[0]) / norm
    return smfn

def smoothFn_windowed_HO(fn, arg, Tf, target=2, a=None, b=None):
    """
    Smoothing function derived for harmonic oscillator. 
    :param fn: Numeric representation of the function to be smoothed.
    :param arg: Functions argument to use for the smoothing.
    :param Tf: Irrelevant as the function is pre-determined.
    :param target: Irrelevant, only taken for compatibility.
    :param a: Irrelevant, only taken for compatibility.
    :param b: Irrelevant, only taken for compatibility.
    :return: Smoothed function values in the same numeric representation.
    """
    # numfile = open("/home/berta/bin/GSTA-heatcap/modules/numHO_300K_10e8_1pk_corr.bin", "rb")
    # sm = pickle.load(numfile)
    # numfile.close()
    A = Tf*const.k*np.square(const.pi)/const.h
    B = 2*Tf*const.k*np.square(const.pi)/const.h
    # Tf = np.float64(300.0)
    smfn = np.empty(shape=len(arg), dtype=np.float64)
    for i in range(len(arg)):
        # W = np.empty(shape=len(arg), dtype=np.float64)
        # W[0:i] = sm[i:0:-1]
        # W[i:] = sm[0:(len(W) - len(sm)) - i]
        W = np.divide(A, np.square(np.cosh(np.multiply(np.subtract(arg, arg[i]), B))))
        norm = integrate.simps(W, dx=arg[1] - arg[0])
        fng = np.multiply(fn, W)
        smfn[i] = integrate.simps(fng, dx=arg[1] - arg[0]) / norm
    return smfn

# daemon for gentle output monitoring
# from daemon.g09daemon import g09daemon

class Calculation():
    """Object for storing the calculation information.

    For each new calculation the program raises an instance of this object. It stores the
    neccessary information for a session and gets updated when its neccessary. The object is
    serialized between runs and for diagnostic purposes as a .CalcFile_ file.
    Attributes:
        freqpath (str): absolute path to the frequency calculation with which the session
        started.
        freqfile (str): the name of the frequency calculation with which the session started.
        workdir (str): absolute path to the directory where freqfile was. Originally the
        serialized object is here.
        mol (Object): cclib object of the molecule read from freqfile.
        Check *http://cclib.github.io/data.html* for further info.
        type (int): specifies the approach to the simulations.
            1. 1 -- MD vibration-wise for chosen number of vibrations.
            2. 2 -- MD with kT/2 energy at all vibrations, chosen number of independent
            trajectories.
        end (int): determines until which point the session should be continued.
            1. 1 -- create velocity generating inputs.
            2. 2 -- run velocity generating inputs.
            3. 3 -- create BOMD inputs.
            4. 4 -- run BOMD inputs.
            5. 5 -- process BOMD data.
        temp (double): the temperature of the calculation.
        target (int): specifies what type of functions will be smoothed.

        velGenInputs (list of str): contains the names of the inputs created generating atomic
        velocities.
        SubCalcs (list of Object): contains the SubCalc objects linked to all investigated
        vibration or trajectory.
        command (str): the command to run Gaussian.
        options (str): options delivered with command after the input.
        maxvib (int): the index of the vibration with maximum frequency still to be
        investigated.
        velGenDir (str): absolute path to the directory where velGenInputs were created.

        NTraj (int) = number of independent trajectories.
        NCalc (int) = number running calculation at a time
        rot (boolean) = option to incorporate rotation in MD simulation
        MDCalcs (list of Object): contains the MDCalc objects linked to independent
        trajectories.
        
    """

    def __init__(self, filePath):
        """Constructor function to Calculation objects.
        Args:
            filePath (str): absolute path to the Gaussian output.
            Originally it should be a frequency calculation output, but it can be any complete
            Gaussian output.
        Initializes the following attribute:
            freqpath = filePath
            freqfile = filePath.split(os.sep)[-1]
            workdir = directory specified in filePath
            mol = processg09output(filePath)
            velGenInputs = []
            SubCalcs = []
            command = None
            options = ""
            maxvib = None
            temp = None
            end = None
            target = 1 # kinetic energy
            type = 1 # vibration-wise approach
            NTraj = None
            rot = False
            VELGENDONE = False
        """
        self.freqpath = filePath
        self.freqfile = filePath.split(os.sep)[-1]
        self.workdir = ""
        for a in filePath.split(os.sep)[:-1]: self.workdir += a + '/'
        self.mol = processg09output(filePath)
        self.moveInputs = []
        self.velGenInputs = []
        self.SubCalcs = []
        self.command = None
        self.options = ""
        self.maxvib = None
        self.temp = None
        self.end = None
        self.target = 1  # kinetic energy
        self.type = None  # vibration-wise approach
        self.NTraj = None
        self.NCalc = 99
        self.rot = False
        self.VELGENDONE = False
        return

    def setCommand(self):
        """Interactive guide to set up the command to run Gaussian. See g09setup for more.
        Not called if given in config file."""
        self.command, self.options = g09setup()

    def setType(self):
        """Interactive guide to set up calculation approach.
        Not called if given in config file."""
        while True:
            try:
                answer = input(
                    "Please choose which approach you wish to follow:\n1 -- perform simulations having energy only on "
                    "one vibration.\n2 -- perform NVE simulations having equal energy on all vibration.\n")
            except KeyboardInterrupt:
                print('\n Exiting GSTA')
                sys.exit(0)
            try:
                self.type = int(answer)
                if self.type == 1 or self.type == 2:
                    break
                else:
                    print("Invalid answer! Options 1 and 2 are available at the moment.")
            except ValueError:
                print("Invalid answer! Options 1 and 2 are available at the moment.")

    def setNTraj(self):
        """Interactive guide to set up the number of independent trajectories.
        Not called if given in config file."""
        while True:
            try:
                answer = input("Please give the number of trajectories you wish to calculate.\n At least 10 "
                               "independent trajectories recommended \n")
            except KeyboardInterrupt:
                print("\n Exiting GSTA")
                sys.exit(0)
            try:
                self.NTraj = int(answer)
                if self.NTraj > 0:
                    break
                else:
                    print("The number of trajectories must be positive!")
            except ValueError:
                print("Invalid answer! Please give a positive integer!")

    def setMaxVib(self):
        """Interactive guide to set up the number of vibrations to be investigated.
        Not called if given in config file."""
        print("The vibrational frequencies are the following (in cm-1):")
        for i in range(len(self.mol.vibfreqs)):
            print(i + 1, " -- ", self.mol.vibfreqs[i])
        print("Please choose a vibrational mode. There will be no MD simmulation for faster modes. 1 -",
              len(self.mol.vibfreqs))
        print("Default is", len(self.mol.vibfreqs))
        maxvib = len(self.mol.vibfreqs)
        while True:
            answer = input("\n")
            if answer == '':
                break
            else:
                try:
                    maxvib = int(answer)
                except ValueError:
                    print("Invalid input. 1 -", len(self.mol.vibfreqs))
                if maxvib not in range(1, len(self.mol.vibfreqs) + 1):
                    print("Invalid number. 1 -", len(self.mol.vibfreqs))
                else:
                    break
        self.maxvib = maxvib

    def setTemp(self):
        """Interactive guide to set up the temperature.
        Not called if given in config file."""
        temp = np.float64(300)
        while True:
            try:
                tempanswer = input("Please give the temperature for the MD simulations in Kelvin.\n") or temp
            except KeyboardInterrupt:
                print("\nExiting GSTA")
                sys.exit(0)
            try:
                temp = np.float64(tempanswer)
                if temp < 0.0:
                    raise ValueError
                else:
                    break
            except ValueError:
                print("Invalid temperature!")
        self.temp = temp

    def setVelGenDir(self, path):
        """Set velGenDir attribute.
        Args:
            path (str): value to be set as velGenDir."""
        self.velGenDir = path

    def addVelGenInput(self, inp):
        """Adds an input name to velGenInputs.
        Args:
            inp (str): name of the input to be added.
        """
        self.velGenInputs.append(inp)

    def addMoveInput(self, inp):
        """Adds an input name to moveInputs.
        Args:
            inp (str): name of the input to be added.
        """
        self.moveInputs.append(inp)

    def setTarget(self):
        """Interactive guide to set up the target type for smoothing.
        Not called if given in config file."""
        # for energies
        # np.float64((np.sqrt(32)*pi*kB*Tf)/h)
        # np.float64((-32*kB**2*Tf**2*pi**3)/(h**2))
        # for coordinates and velocities
        # np.float64((np.sqrt(8)*pi*kB*Tf)/h)
        # np.float64((-8*kB**2*Tf**2*pi**3)/(h**2))
        while True:
            choose = 0
            choose = input(
                "\nPlease choose the constants A and B according to the function you wish to smooth:\n1. Type '1' for "
                "constants used for filtering energies (this is the default).\n2. Type '2' for constants used for "
                "filtering coordinates or velocities.\n3. Type '3' for giving the constants manually.\n4. Type '?' "
                "for information about the pre-defined constants.\n") or 1
            if choose == 1:
                break
            elif choose == 2:
                self.target = 2
                break
            elif choose == 3:
                print("The function used is given by the formula:\ng(tau) = A * exp(B * tau**2)\n")
                a = 1
                a = input("Please give A in numpy format!\nnumpy imported as np\n")
                try:
                    self.A = np.float64(a)
                except ValueError:
                    print("Input cannot be interpreted.")
                b = 1
                b = input("Please give B in numpy format!\nnumpy mported as np\n")
                try:
                    self.B = np.float64(b)
                    break
                except ValueError:
                    print("Input cannot be interpreted.")
            elif choose == '?':
                print("The function used is given by the formula:\ng(tau) = A * exp(B * tau**2)\n")
                print("By default, for energies:\n A = (sqrt(32)*pi*kB*Tf)/h\n B = (-32*kB**2*Tf**2*pi**3)/(h**2)\n")
                print(
                    "For coordinates and velocities:\n A = (sqrt(8)*pi*kB*Tf)/h\n B = (-8*kB**2*Tf**2*pi**3)/(h**2)\n")
            else:
                print("Please choose from the given options!\n")

    def setEnd(self):
        """Interactive guide to set up the end of the session.
        Not called if given in config file."""
        print("Please choose until which point you wish to proceed automatically.")
        while True:
            answer = 0
            print(
                "1 -- create velocity generating inputs.\n2 -- run velocity generating inputs.\n3 -- create BOMD "
                "inputs.\n4 -- run BOMD inputs.\nProcessing has to be done after simulations are done, by calling "
                "'GSTA-hc -t pp -f .CalcFile*'")
            try:
                answer = int(input("\n"))
            except KeyboardInterrupt:
                print("\n Exiting GSTA")
                sys.exit(0)
            except ValueError:
                pass
            if answer in range(1, 6):
                self.end = answer
                break
            else:
                print("Please choose from the given options!\n")
        # if answer == 5: self.setTarget() # not stable, taken out
        return

    def setAllInfo(self):
        """Interactive guide to set up command, maxvib, temp and end attributes."""
        self.setCommand()
        self.setMaxVib()
        self.setTemp()
        self.setEnd()

    def addSubCalc(self, inp):
        """Creates a SubCalc object for one vibration or trajectory and adds it to the SubCalcs
        list. Substitute for addVibCalc() and addMDCalc().
        See SubCalc.__init__() for more.
        Args:
            inp (str): name of the input file.
        """
        self.SubCalcs.append(SubCalc(parent=self, inp=inp))

    def updateSubCalc(self, out):
        """Updates the corresponding SubCalc object after a velgen calculation has terminated.
        Args:
            out (str): absolute path to the terminated calculation.

        Velocities are added to the mol object.
        """
        if "_velgen." in out:
            name = "traj" + out.split("_traj")[1].split('_velgen.')[0]
        else:
            name = "traj" + out.split("_traj")[1].split('.')[0]
        for SC in self.SubCalcs:
            if SC.name == name:
                SC.mol = processg09output(out, isBOMD=True)
                SC.mol.atomcoords = 0.52917724900001 * getCoords(out)[0]
                if "_velgen." in out:
                    SC.VELGENDONE = True
                elif SC.rot:
                    SC.VELGENDONE = True
                else:
                    SC.MOVED = True
                    SC.mol.modeEkin = []
                    with open(out, 'r') as f:
                        READ = False
                        for line in f:
                            if "Summary of normal mode sampling:" in line:
                                READ = False
                                break
                            if READ and len(line.split()) == 4:
                                if line.split()[0] == "Mode":
                                    SC.mol.modeEkin.append((np.float64(
                                        line.split()[3].split('D')[0] + 'E' + line.split()[3].split('D')[
                                            1]) ** 2) * 3.34642001457887E-28)
                            if "MW displacement        MW velocity" in line:
                                READ = True
                if SC.VELGENDONE:
                    SC.workdir = self.workdir + os.sep + "MD"
                logging.info("PID" + str(os.getpid()) + ": " + SC.name + " is updated")
                break
        return

    def removeSubCalc(self, out):
        """Removes the corresponding SubCalc object after a velgen calculation has failed.
                Args:
                    out (str): absolute path to the terminated calculation.
        """
        if self.type == 1:
            name = "vib" + out.split("_vib")[1].split('.')[0]
        elif self.type == 2:
            name = "traj" + out.split("_traj")[1].split('.')[0]
        for SC in self.SubCalcs:
            if SC.name == name:
                self.SubCalcs.remove(SC)
                logging.info("PID" + str(os.getpid()) + ": " + SC.name + " is removed")
                break
        return

    def isVELGENDONE(self):
        """ Loops over the subcalculations and returns True if all velocity generating jobs have been done.
        Returns (boolean): True if all calculations are done, False otherwise.
        """
        for SC in self.SubCalcs:
            if not SC.VELGENDONE:
                # print(SC.name,SC.VELGENDONE)
                return False
        logging.info("PID" + str(os.getpid()) + ": " + str(
            len(self.SubCalcs)) + " velgen calculations are done, proceeding to MD.")
        return True

    def isMOVED(self):
        """ Loops over the subcalculations and returns True if all coordinate transformation jobs have been done.
        Returns (boolean): True if all calculations are done, False otherwise.
        """
        for SC in self.SubCalcs:
            if not SC.MOVED:
                return False
        logging.info("PID" + str(os.getpid()) + ": " + str(
            len(self.SubCalcs)) + " geometry transform calculations are done, proceeding to velgen.")
        return True

    def save(self, lock=None):
        """Serializes itself to `".CalcFile_"+self.freqfile.split('.')[0]`
        Args:
            lock (Object): a filelock object can be given to release locked files after changes.
            Otherwise the file is locked to avoid multiple reading and/or writing
            attepmts.
        """
        os.chdir(self.workdir)
        while True:
            try:
                if lock == None:
                    lock = fl.FileLock(".CalcFile_" + self.freqfile.split('.')[0])
                    lock.acquire()
                cf = open(".CalcFile_" + self.freqfile.split('.')[0], 'wb')
                pickle.dump(self, cf)
                logging.info("PID" + str(os.getpid()) + ": Calc is written out to " + cf.name)
                cf.close()
                break
            except fl.FileLockException:
                logging.warn("PID" + str(os.getpid()) + ": " + cf.name + " lock cannot be acquired.")
                sys.exit(0)
            finally:
                lock.release()
        return

    def load(self, expectedObjPath):
        """Loads the serialized object.

        Args:
            expectedObjPath (str): absolute path to the serialized object.
        Returns:
            Object: the loaded object.
            Object: the filelock Object to the checkpoint file.
        Raises:
            IOError: if the file cannot be load.

        The function locks the accessed file and returns the lock. After changes taken place,
        lock can be released by saving the file by save(lock) or by releasing it by
        lock.release()).
        """
        if os.path.isfile(expectedObjPath):
            while True:
                try:
                    tree = expectedObjPath.split(os.sep)
                    d = ""
                    for i in range(len(tree) - 1):
                        d += os.sep + tree[i]
                    os.chdir(d)
                    lock = fl.FileLock(expectedObjPath.split(os.sep)[-1])
                    lock.acquire()
                    CO = open(expectedObjPath, 'rb')
                    LO = pickle.load(CO)
                    CO.close()
                    break
                except fl.FileLockException:
                    logging.warn("PID" + str(os.getpid()) + ": " + CO.name + " lock cannot be acquired.")
                    sys.exit(0)
        else:
            logging.warn("Pickled Calc File not found at " + expectedObjPath)
            logging.warn("Exiting")
            raise IOError

        return LO, lock

    def updateMDneg(self, out):
        """Updates the corresponding SubCalc object after an MD simulation has terminated.
        Supposed to get a calculation started with negative velocities which is implied by
        '_neg'.
        Args:
            out (str): absolute path to the terminated calculation.

        MDneg = True for the corresponding SubCalc and MDoutput2 is set.
        """
        os.chdir(self.workdir)
        if self.type == 1:
            name = "vib" + out.split("_vib")[1].split('.')[0].split("_MD_")[0]
        elif self.type == 2:
            name = "traj" + out.split("_traj")[1].split('.')[0].split("_MD_")[0]
        for SC in self.SubCalcs:
            if SC.name == name:
                SC.MDneg = True
                SC.MDoutput2 = out
                logging.info("PID" + str(os.getpid()) + ": " + SC.name + " updated: MDneg = " + str(SC.MDneg))
                break
        return

    def updateMDpos(self, out):
        """Updates the corresponding SubCalc object after an MD simulation has terminated.
        Supposed to get a calculation started with original velocities which is implied by
        '_pos'.
        Args:
            out (str): absolute path to the terminated calculation.

        MDpos = True for the corresponding SubCalc and MDoutput1 is set.
        """
        os.chdir(self.workdir)
        if self.type == 1:
            name = "vib" + out.split("_vib")[1].split('.')[0].split("_MD_")[0]
        elif self.type == 2:
            name = "traj" + out.split("_traj")[1].split('.')[0].split("_MD_")[0]
        for SC in self.SubCalcs:
            if SC.name == name:
                SC.MDpos = True
                SC.MDoutput1 = out
                logging.info("PID" + str(os.getpid()) + ": " + SC.name + " updated: MDpos = " + str(SC.MDpos))
                break
        return

    def processConfig(self, configfile):
        if os.path.isfile(configfile):
            with open(configfile, 'r') as config:
                for line in config:
                    if "command" in line.split('=')[0]:
                        self.command = line.split('=')[1].strip()
                    elif "maxvib" in line.split('=')[0]:
                        if line.split('=')[1].strip() == "all":
                            self.maxvib = len(self.mol.vibfreqs)
                        else:
                            self.maxvib = int(line.split('=')[1].strip())
                    elif "temp" in line.split('=')[0]:
                        self.temp = np.float64(line.split('=')[1].strip())
                    elif "end" in line.split('=')[0]:
                        self.end = int(line.split('=')[1].strip())
                    elif "target" in line.split('=')[0]:
                        self.target = int(line.split('=')[1].strip())
                    elif "type" in line.split('=')[0]:
                        self.type = int(line.split('=')[1].strip())
                    elif "NTraj" in line.split('=')[0]:
                        self.NTraj = int(line.split('=')[1].strip())
                    elif "options" in line.split('=')[0]:
                        self.options = line.split('=')[1].strip()
                    elif "NCalc" in line.split('=')[0]:
                        self.NCalc = int(line.split('=')[1].strip())
                    elif "rotation" in line.split('=')[0]:
                        if line.split('=')[1].strip() == "true":
                            self.rot = True
                        elif line.split('=')[1].strip() == "false":
                            self.rot = False
        else:
            logging.warn("PID" + str(os.getpid()) + ": config file cannot be accessed.")
        return


class SubCalc(Calculation):
    """Object for storing information about one thread of the session. It is a child of a
    Calculation object, and stored a SubCalcs list of one.
    Attributes:
        name (str): vibX where X is the index of vibration, or trayX where X is the index of the
        independent trajectory.
        freqfile (str): original freq output. See Calculation docstring.
        type (int): type of parent. 1 is for vibvise sampling, 2 is for NVE.
        workdir (str): subdirectory of parent, velgen or MD, depending which state the session
        is in.
        mol (Object): mol of parent. See Calculation docstring.
        command (str): command of parent. See Calculation docstring.
        options (str): options of parent. See Calculation docstring.
        temp (double): temp of parent. See Calculation docstring.
        end (int): end of parent. See Calculation docstring.
        target (int): target of parent. Specifies the smoothing function.
        inputname (str): corresponding velGenInput.
        MDinput1 (str): input name for MD simulation ending '_pos'.
        MDinput2 (str): input name for MD simulation ending '_neg'.
        MDoutput1 (str): output path to MD simulation ending '_pos', if it is done.
        MDoutput2 (str): output path to MD simulation ending '_neg', if it is done.
        vibFreq (double): vibrational frequency. The smallest if sampling type is NVE.
        MP (int): MaxPoints option for gaussian. 200+decay of the smoothing function+1
        VELGENDONE (bool): indicates if velgen calculation is done.
        BOMDDONE (bool): indicates if both BOMD calculations are done.
        MDpos (bool): indicates if BOMD calculation ending '_pos' is done.
        MDneg (bool): indicates if BOMD calculation ending '_neg' is done.
      
        vibIndex (int): index of vibrational frequency (starts with 0). Only for vibwise
        calculation.

        rand (int): Random seed for NVE sampling.
    """

    def __init__(self, parent, inp):
        """Constructor function to SubCalc objects.
        Args:
            parent (Object): Calculation object of the session.
            inp (str): name of the input file created generating velocities.
        Initializes the following attribute:
            name = trajX
            freqfile = parent.freqfile
            workdir = parent.velGenDir
            mol = parent.mol
            command = parent.command
            options = parent.options
            temp = parent.temp
            end = parent.end
            rot = parent.rot
            inputname = inp
            VELGENDONE = False
            BOMDDONE = False
            MDpos = False
            MDneg = False
            only for NVE sampling:
            rand = None
            MDinput1 = MDinput2 = None
        """
        self.freqfile = parent.freqfile
        self.workdir = parent.velGenDir
        self.mol = parent.mol
        self.command = parent.command
        self.options = parent.options
        self.temp = parent.temp
        self.end = parent.end
        self.inputname = inp
        self.VELGENDONE = False
        self.BOMDDONE = False
        self.MDpos = False
        self.MDneg = False
        self.name = "traj" + inp.split("_traj")[1].split('.')[0]
        self.velgenInput = None
        self.rand = None
        self.MDinput1 = self.MDinput2 = None
        self.rot = parent.rot
        if not parent.rot:
            self.MOVED = False

    def save(self, lock=None):
        """Serializes itself to `".SubCalcFile_"+self.name`
        Args:
            lock (Object): a filelock object can be given to release locked files after changes.
            Otherwise the file is locked to avoid multiple reading and/or writing
            attempts.
        """
        os.chdir(self.workdir)
        while True:
            try:
                if lock == None:
                    lock = fl.FileLock(".SubCalcFile_" + self.name)
                    lock.acquire()
                cf = open(".SubCalcFile_" + self.name, 'wb')
                pickle.dump(self, cf)
                logging.info("PID" + str(os.getpid()) + ": SubCalc " + self.name + " is written out to " + cf.name)
                cf.close()
                break
            except fl.FileLockException:
                sleep(random.uniform(1, 2))
            finally:
                lock.release()
        return


######################################################################
#                               TESTS                                #
######################################################################
def main():
    dirpath = "/media/berta/kaptar_home/calculations/heat_capacity/MeOH/MD/"
    # filepath1 = dirpath + "MeOH_traj1_MD_pos.out"
    # filepath2 = dirpath + "MeOH_traj1_MD_neg.out"

    # mol = processg09output(filepath, isBOMD=True)
    # createMDInput(mol,filepath)
    # createVelInputs(mol,filepath)

    # trajecs1 = getTrajSummary(filepath1)
    # trajecs2 = getTrajSummary(filepath2)

    # xyz = getCoords(filepath)
    # t_pos = processTrajSummary(trajecs1)
    # t_neg = processTrajSummary(trajecs2)
    # t = np.multiply(np.concatenate((-t_neg[0][:0:-1], t_pos[0])), 1e-15)

    # Ep_pos = getEpot(filepath1)
    # Ep_neg = getEpot(filepath2)
    # Ep = np.concatenate((Ep_neg[0][:0:-1], Ep_pos[0]))

    t = np.linspace(-5e-13, 5e-13, 1e4+1)
    # nu_list = np.linspace(1.2e13, 1.2e14, 10)
    # T_list = np.linspace(100.0, 1000.0, 10)
    # for nu in nu_list:
    #     for T in T_list:
    T = 300.0
    nu2 = 1.2e13
    nu1 = 5e12
    phi = 2e13
    # v = np.cos(t*2*const.pi*nu1)+np.cos(t*2*const.pi*nu2+phi)
    v = np.cos(t * 2 * const.pi * nu1)

    # smv = smoothFn_windowed(Ep, t, 300)
    smv2 = smoothFn_windowed_HO(v, t, T)
    plt.plot(t, v, "o", t, smv2, "x")
    plt.show()

    ac = autocorr_manual(v)
    smac = autocorr_manual(smv2)
    # plt.plot(t, ac)
    # plt.show()

    dos = getDos(ac)
    smdos = getDos(smac)
    dosQM = getDosQM(dos, T)
    # f, axarr = plt.subplots(2)
    # axarr[0].plot(range(6000), dos)
    # axarr[1].plot(range(6000), dosQM)
    # axarr[0].set_title("dos")
    # axarr[1].set_title("dosQM")
    # plt.plot(range(6000), dos, range(6000), dosQM, range(6000), smdos)
    # plt.show()
    with open("/home/berta/num-dos", "w") as f:
        for i in range(6000):
            f.write(str(i)+" "+str(dos[i])+" "+str(dosQM[i])+"\n")

    int_class = integrate.simps(dos, dx=1.0)
    int_QM = integrate.simps(dosQM, dx=1.0)
    int_sm = integrate.simps(smdos, dx=1.0)
    cv_1PT = 1 * const.R * int_QM / int_class
    cv_sm = 1 * const.R * int_sm / int_class
    cv_exact1 = (const.R * np.square((const.h * nu1) / (2 * T * const.k))) / (np.square(np.sinh((const.h * nu1) / (2 * T * const.k))))
    # cv_exact2 = (const.R * np.square((const.h * nu2) / (2 * T * const.k))) / (np.square(np.sinh((const.h * nu2) / (2 * T * const.k))))
    cv_exact = cv_exact1 # + cv_exact2
    E = np.square(v)
    Esm = np.square(smv2)
    E_avr = (1 / ((len(E) - 1) * 1e16)) * integrate.simps(E, dx=1e16)
    Esm_avr = (1 / ((len(Esm) - 1) * 1e16)) * integrate.simps(Esm, dx=1e16)
    cv_sm2 = 1 * const.R * Esm_avr / E_avr

    print(cv_exact, cv_1PT, cv_sm, cv_sm2)
    # with open("/home/denes/numtest","a") as f:
    #     f.write("Numeric test 1 ps, 0.1 fs timestep, nu: "+str(nu)+" T: "+str(T)+"\n")
    #     f.write("autocorrelation function\n")
    #     for n in ac: f.write(str(n))
    #     f.write("\nsmoothed autocorrelation function\n")
    #     for n in smac: f.write(str(n))
    #     f.write("\nDoS\n")
    #     for n in dos: f.write(str(n))
    #     f.write("\nDoS weighted\n")
    #     for n in dosQM: f.write(str(n))
    #     f.write("\nDoS from smoothed velocity\n")
    #     for n in smdos: f.write(str(n))
    #     f.write("\n")
    # with open("/home/denes/numtest_cv","a") as f:
    #     f.write("Numeric test 1 ps, 0.1 fs timestep, nu: " + str(nu) + " T: " + str(T) + "\n")
    #     f.write("QHO: "+str(cv_exact)+" Smoothing DoS: "+str(cv_sm)+" Smoothing: "+str(const.R*Esm_avr/E_avr)+" 1PT: "+str(cv_1PT)+"\n")


if __name__ == "__main__":
    main()
