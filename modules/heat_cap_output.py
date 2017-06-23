#!/usr/bin/env python3
import numpy as np
import modules.g09BOMD_filter
from daemon.g09daemon import tail
import modules.Heatcapacity_subproc as hcs
import modules.GetHeatcapacity_util as hcu
import modules.g09BOMD_filter as g09f
# import multiprocessing.pool
from multiprocessing import Pool
from functools import partial
import glob, sys, os, pickle, code

try:
    import matplotlib.pyplot as plt
except ImportError:
    pass
import scipy.constants as const
import scipy.integrate as integrate


# class NoDaemonProcess(multiprocessing.Process):
#     # make 'daemon' attribute always return False
#     def _get_daemon(self):
#         return False
#     def _set_daemon(self, value):
#         pass
#     daemon = property(_get_daemon, _set_daemon)
#
# class Pool(multiprocessing.pool.Pool):
#     Process = NoDaemonProcess

class CVoutputException(Exception):
    pass


class CV_output():
    """
    Object to handle processed trajectories and print the output.
    """

    def __init__(self, Calc):
        self.workdir = Calc.workdir
        self.SCfiles = glob.glob(os.path.join(self.workdir, "MD", ".SubCalcFile*"))
        self.cvout = os.path.join(self.workdir, Calc.freqfile + ".cv.out")
        with open(self.cvout,'w') as f:
            self.printHeader(f)

    #
    # def __del__(self):
    #     self.cvout.close()

    def printHeader(self,FH):
        """ Writes the basic software information into the output file."""
        header = "Heat capacity calculation via trajectory smoothing. Born-Oppenheiemer molecular dynamics performed as implemented into Gaussian software. Check https://www.gaussian.com for details." + os.linesep + os.linesep + "Data process is done as proposed by Madarasz et al. in **insert publication here**" + os.linesep + os.linesep + "The software employs python3 and cclib utilities and the source code is available at https://github.com/bertadenes/g09entropy.git" + os.linesep + os.linesep + "Developed by" + os.linesep + "D. Berta    berta.denes@ttk.mta.hu*" + os.linesep + "D. Ferenc   ferenc.david@ttk.mta.hu" + os.linesep + "A. Madarasz madarasz.adam@ttk.mta.hu" + os.linesep + "Reported issues regarding the software are highly appreciated." + os.linesep + os.linesep + "The developers acknowledge the help of Vilmos Nagy in programing issues."
        try:
            FH.write(header)
            FH.write(os.linesep)
        except:
            raise CVoutputException
        return

    def print(self,line):
        with open(self.cvout,'a') as f:
            f.write(line)
            f.write(os.linesep)
        return

class CV_postProcess():
    """ Utility to process a set of MD calculation carried out using entropy. Please note that only applicable for
        type = 2 calculations at the moment.
    """

    def __init__(self, Calc, workdir):
        self.workdir = workdir
        self.temp = Calc.temp
        try:
            if Calc.rot:
                self.df = (3 * len(Calc.mol.atomnos) - 3)
            else:
                self.df = (3 * len(Calc.mol.atomnos) - 6)
        except AttributeError:
            self.df = (3 * len(Calc.mol.atomnos) - 6)
        self.trajs = []
        self.wrongtrajs = []
        self.rndtemps = []
        self.rtemps = []
        self.outputsPOS = []
        self.outputsNEG = []
        self.dos = []
        self.dosQM = []
        self.dosG = []
        self.cv = []
        self.cvG = []
        self.Ekinsmavr = []
        self.Epotsmavr = []
        self.Tclass = []
        for SC in Calc.SubCalcs:
            out1ok = out2ok = False
            if SC.BOMDDONE:
                out1ok = out2ok = True
            elif SC.MDinput1 is not None and SC.MDinput2 is not None:
                if len(glob.glob(os.path.join(workdir, SC.MDinput1.split('.')[0] + "*log"))) == 1:
                    SC.MDoutput1 = glob.glob(os.path.join(workdir, SC.MDinput1.split('.')[0] + "*log"))[0]
                elif len(glob.glob(os.path.join(workdir, SC.MDinput1.split('.')[0] + "*out"))) == 1:
                    SC.MDoutput1 = glob.glob(os.path.join(workdir, SC.MDinput1.split('.')[0] + "*out"))[0]
                else:
                    print(SC.name + ": there is no output for MDinput1 = " + SC.MDinput1)
                    continue
                if len(glob.glob(os.path.join(workdir, SC.MDinput2.split('.')[0] + "*log"))) == 1:
                    SC.MDoutput2 = glob.glob(os.path.join(workdir, SC.MDinput2.split('.')[0] + "*log"))[0]
                elif len(glob.glob(os.path.join(workdir, SC.MDinput2.split('.')[0] + "*out"))) == 1:
                    SC.MDoutput2 = glob.glob(os.path.join(workdir, SC.MDinput2.split('.')[0] + "*out"))[0]
                else:
                    print(SC.name + ": there is no output for MDinput2 = " + SC.MDinput2)
                    continue
            else:
                print(
                    SC.name + ": either one or two MD inputs are not assigned to this trajectory. Check if they have been created.")
            try:
                with open(SC.MDoutput1, 'rb') as outfile:
                    if "termination" in tail(outfile):
                        if "Error" in tail(outfile):
                            print(SC.MDoutput1 + " has terminated with an error.")
                        else:
                            out1ok = True
                with open(SC.MDoutput2, 'rb') as outfile:
                    if "termination" in tail(outfile):
                        if "Error" in tail(outfile):
                            print(SC.MDoutput2 + " has terminated with an error.")
                        else:
                            out2ok = True
            except AttributeError:
                pass
            if out1ok and out2ok:
                self.trajs.append(SC.name)
                self.rndtemps.append(SC.rndtemp)
                try:
                    self.rtemps.append(SC.rtemp)
                except:
                    self.rtemps.append(0)
                self.outputsPOS.append(SC.MDoutput1)
                self.outputsNEG.append(SC.MDoutput2)
            else:
                self.wrongtrajs.append(SC.name)

    def removeWrongTrajs(self):
        for t in self.wrongtrajs:
            try:
                index = self.trajs.index(t)
                self.trajs.remove(t)
                self.rndtemps.remove(self.rndtemps[index])
                self.outputsNEG.remove(self.outputsNEG[index])
                self.outputsPOS.remove(self.outputsPOS[index])
            except ValueError:
                pass

    def getCV_fromDos(self, traj, vstep=1.0):
        """Integrates the DoS over the frequency domain reslting in the heat capacity in J/K*mol"""
        index = self.trajs.index(traj)
        int_class = integrate.simps(self.dos[index], dx=vstep)
        int_QM = integrate.simps(self.dosQM[index], dx=vstep)
        int_G = integrate.simps(self.dosG[index], dx=vstep)
        self.cv.append(self.df * const.R * int_QM / int_class)
        self.cvG.append(self.df * const.R * int_G / int_class)
        # print(self.df,int_QM,int_class,self.df * const.R * int_QM / int_class)
        return

    def getCV_fromVelSmoothing(self):
        return hcu.fitLinear(self.Tclass,np.multiply(np.array(self.Ekinsmavr),2))

    def getCV_fromTrajSmoothing(self):
        return hcu.fitLinear(self.Tclass,np.add(np.array(self.Ekinsmavr),np.array(self.Epotsmavr)))

def calcDOS(traj, ppobj):
    """Not function of CV_postProcess because it raises error when called by Pool."""
    index = ppobj.trajs.index(traj)
    try:
        dos = hcs.DOS_from_velocity(ppobj.outputsPOS[index], ppobj.outputsNEG[index], plot=False, temp=ppobj.temp)
        ppobj.dos.append(dos)
        # print(traj+" added. DOS length:",len(ppobj.dos))
    except g09f.NonSmoothableTraj:
        ppobj.wrongtrajs.append(traj)
    return

def calcSmoothedEnergies(traj, ppobj):
    """Not function of CV_postProcess because it raises error when called by Pool."""
    index = ppobj.trajs.index(traj)
    try:
        Ep, Ek, T = hcs.get_Ekinsmavrg_Epotsmavrg(ppobj.outputsPOS[index], ppobj.outputsNEG[index], temp=ppobj.temp)
        ppobj.Ekinsmavr.append(Ek)
        ppobj.Epotsmavr.append(Ep)
        ppobj.Tclass.append(T)
    except g09f.NonSmoothableTraj:
        ppobj.wrongtrajs.append(traj)
    return

def mp_calcDOS(ppobj):
    """creates daemon that cannot use this function due it'd create daemons too"""
    with Pool(processes=8) as pool:
        pool.map(partial(calcDOS, ppobj=ppobj), ppobj.trajs)
    return

def postProcess(CF_name):
    CF = open(CF_name, "rb")
    Calc = pickle.load(CF)
    workdir = os.path.join(os.path.dirname(os.path.realpath(CF_name)), "MD")
    CF.close()
    cvout = CV_output(Calc)
    if Calc.type == 2:
        pp = CV_postProcess(Calc, workdir)
        for traj in pp.trajs:
            calcDOS(traj, pp)
            pp.dosQM.append(hcu.getDosQM(pp.dos[-1], pp.temp))
            pp.dosG.append(hcu.getDosG(pp.dos[-1], pp.temp))
            calcSmoothedEnergies(traj, pp)
            cvout.print("{0:s} processed".format(traj))
            # break
        pp.removeWrongTrajs()
        for traj in pp.trajs:
            pp.getCV_fromDos(traj)
        CV_v, err_v = pp.getCV_fromVelSmoothing()
        CV_t, err_t = pp.getCV_fromTrajSmoothing()
        #     break
        # f, axarr = plt.subplots(3)
        # axarr[0].plot(range(6000), pp.dos[0])
        # axarr[1].plot(range(6000), pp.dosG[0])
        # axarr[2].plot(range(6000), pp.dosQM[0])
        # axarr[0].set_title("dos")
        # axarr[1].set_title("dosG")
        # axarr[2].set_title("dosQM")
        # plt.show()
        cvout.print("Trajectory information and smoothed energies")
        cvout.print("traj\tvib temp\trot temp\tkin temp\tmean Ekinsm\tmean Epotsm A.U.")
        for i in range(len(pp.trajs)):
            cvout.print("{0:s}\t{1:f}\t{2:d}\t{3:f}\t{4:f}\t{5:f}".format(pp.trajs[i], pp.rndtemps[i], pp.rtemps[i], pp.Tclass[i], pp.Ekinsmavr[i] / 627503, pp.Epotsmavr[i] / 627503.0))
        cvout.print("Heat capacities:")
        cvout.print("Berens: {0:f} +- {1:f} cal/Kmol".format(np.mean(pp.cv)/4.184, np.std(pp.cv)/4.184))
        cvout.print("Gauss: {0:f} +- {1:f} cal/Kmol".format(np.mean(pp.cvG)/4.184, np.std(pp.cvG)/4.184))
        cvout.print("Velocity smoothing: {0:f} +- {1:f} cal/Kmol".format(CV_v, err_v,))
        cvout.print("Trajectory smoothing: {0:f} +- {1:f} cal/Kmol".format(CV_t, err_t))
        cvout.print("Based on {0:d} trajectories".format(len(pp.trajs)))
        cvout.print("{0:d} wrong trajectories".format(len(pp.wrongtrajs)))
        cvout.print("{0:d} Density of states calculated".format(len(pp.dos)))


def main():
    postProcess(".CalcFile_H2O_B3PW91")

    # dos = hcs.DOS_from_velocity("MD/water_traj1_MD_pos.out","MD/water_traj1_MD_neg.out",plot=False,temp=100)
    # cv = hcu.getCV_fromDos(dos)
    # print("classical cv",cv,1.9872036*6)
    # cvsm = hcu.getCV_fromDos(dossm)
    # dosQM = hcu.getDosQM(dos,700.0)
    # dosG = hcu.getDosG(dos,700.0)
    # cvQM = hcu.getCV_fromDos(dosQM)
    # cvG = hcu.getCV_fromDos(dosG)
    # print("QM weighted cv",cvQM,1.9872036*6*cvQM/cv)
    # print("Gauss weighted cv",cvG, 1.9872036 * 6 * cvG / cv)
    # print("Gauss smoothed cv",cvsm, 1.9872036 * 6 * cvsm / cv)
    # f, axarr = plt.subplots(2, 2)
    # axarr[0, 0].plot(range(6000), dos)
    # axarr[1, 0].plot(range(6000), dossm)
    # axarr[1, 1].plot(range(6000), dosG)
    # axarr[0, 1].plot(range(6000), dosQM)
    # axarr[0, 0].set_title("dos")
    # axarr[1, 0].set_title("dossm")
    # axarr[1, 1].set_title("dosG")
    # axarr[0, 1].set_title("dosQM")
    # plt.show()


if __name__ == "__main__":
    main()
