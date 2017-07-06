#!/usr/bin/env python3
"""
Utility to perform heat capacity calculations using g09 Born-Oppenheimer Molecular Dynamics and smoothing the data with Gaussian functions.

Developed by
D. Berta    berta.denes@ttk.mta.hu
D. Ferenc   ferenc.david@ttk.mta.hu
A. Madarasz madarasz.adam@ttk.mta.hu

Dependencies
python3 is recommended and some additional associated package:
python3-cclib may require openbabel
python3-numpy
python3-scipy
"""

programName = "GSTA-hc"

from daemon.g09daemon import g09daemon
import modules.g09BOMD_filter as g
import sys, os, pickle, logging, time, argparse, glob
import numpy as np
import modules.Heatcapacity_subproc as hcs
import modules.heat_cap_output as hcout


def main():
    BOMDDONE = VELGENDONE = VELGENWRONG = MOVE = False
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file')
    parser.add_argument('-c', '--config')
    parser.add_argument('-t', '--type')
    args = parser.parse_args()

    if args.type == "pp":
        if args.file != None:
            if os.fork():
                sys.exit(0)
            hcout.postProcess(args.file)
            sys.exit(0)
        else:
            hcout.main()
            sys.exit(0)

    if args.file != None:
        filepath = args.file
    else:
        filepath = None
        print("Please give the path to the output!")
        try:
            filepath = input("\n") or filepath
            # if filepath is None:
            #     print("Please give the path to the output!")
            #     filepath = input("\n")
        except KeyboardInterrupt:
            print("\n Exiting GSTA")
            sys.exit(0)
        except TypeError:
            # print("Please give the path to the output!")
            pass
    if os.path.isfile(filepath):
        if not os.path.isabs(filepath):
            filepath = os.path.realpath(filepath)
    else:
        print("File not found. " + filepath)
        sys.exit(0)

    Calc = g.Calculation(filepath)
    if args.type == "move":
        MOVE = True
    elif args.type == "velgen":
        VELGENDONE = True

    elif args.type == "MD":
        BOMDDONE = True

    elif args.type == "velgenwrong":
        VELGENDONE = True
        VELGENWRONG = True

    if BOMDDONE:
        fullpath = os.path.realpath(args.file)
        tree = fullpath.split(os.sep)
        expectedObjPath = ""
        for i in range(len(tree) - 2): expectedObjPath += tree[i] + os.sep
        chks = glob.glob(expectedObjPath + os.sep + ".CalcFile_*")
        if len(chks) == 0:
            logging.warn("Pickled Calc File not found at " + expectedObjPath)
            logging.warn("Exiting")
            sys.exit(0)
        elif len(chks) == 1:
            expectedObjPath = chks[0]
        else:
            logging.warn(
                "More then one pickled Calc File found at " + expectedObjPath + ". Please remove those not linked to " + args.file)
            logging.warn("Exiting")
            sys.exit(0)
        try:
            Calc, FL = Calc.load(expectedObjPath)
        except IOError:
            sys.exit(0)

        os.chdir(Calc.workdir)
        logging.basicConfig(filename=Calc.freqfile.split('.')[0] + ".lg", level=logging.INFO,
                            format='%(asctime)s %(message)s')
        logging.info("PID" + str(os.getpid()) + ": " + programName + " started with VELGENDONE = " + str(
            VELGENDONE) + " and BOMDDONE = " + str(BOMDDONE))
        logging.info("PID" + str(os.getpid()) + ": Calc is load from " + expectedObjPath)

        if "_MD_neg" in filepath.split(os.sep)[-1]:
            Calc.updateMDneg(filepath)

        if "_MD_pos" in filepath.split(os.sep)[-1]:
            Calc.updateMDpos(filepath)

        name = "traj" + filepath.split("_traj")[1].split('.')[0].split("_MD_")[0]
        for SC in Calc.SubCalcs:
            if SC.name == name:
                if SC.MDpos and SC.MDneg:
                    SC.BOMDDONE = True
                    logging.info("PID" + str(os.getpid()) + ": " + SC.name + " is ready for processing.")
                    if SC.type == 1:
                        # vibwise process
                        SC.save()
                        FL.release()
                        hcs.heatcapacity_processing(SC)
                    elif SC.type == 2:
                        # nve process
                        FL.release()

                elif not SC.MDpos and SC.MDneg:
                    logging.info("PID" + str(os.getpid()) + ": " + SC.name + '.' + "MDpos = " + str(
                        SC.MDpos) + ". Waiting for it, ending thread.")
                    FL.release()
                    sys.exit(0)
                elif SC.MDpos and not SC.MDneg:
                    logging.info("PID" + str(os.getpid()) + ": " + SC.name + '.' + "MDneg = " + str(
                        SC.MDneg) + ". Waiting for it, ending thread.")
                    FL.release()
                    sys.exit(0)
                else:
                    logging.warn("PID" + str(
                        os.getpid()) + ": both " + SC.name + '.' + "MDpos and " + SC.name + '.' + "MDneg are False. Something went wrong. Are the filenames correct?")
                    logging.info("PID" + str(os.getpid()) + ": Exiting")
                    FL.release()
                    sys.exit(0)

    elif VELGENDONE:
        fullpath = os.path.realpath(args.file)
        tree = fullpath.split(os.sep)
        expectedObjPath = ""
        for i in range(len(tree) - 2): expectedObjPath += tree[i] + os.sep
        chks = glob.glob(expectedObjPath + os.sep + ".CalcFile_*")
        if len(chks) == 0:
            logging.warn("Pickled Calc File not found at " + expectedObjPath)
            logging.warn("Exiting")
            sys.exit(0)
        elif len(chks) == 1:
            expectedObjPath = chks[0]
        else:
            logging.warn(
                "More then one pickled Calc File found at " + expectedObjPath + ". Please remove those not linked to " + args.file)
            logging.warn("Exiting")
            sys.exit(0)
        try:
            Calc, FL = Calc.load(expectedObjPath)
        except IOError:
            sys.exit(0)

        os.chdir(Calc.workdir)
        logging.basicConfig(filename=Calc.freqfile.split('.')[0] + ".lg", level=logging.INFO,
                            format='%(asctime)s %(message)s')
        logging.info("PID" + str(os.getpid()) + ": " + programName + " started with VELGENDONE = " + str(
            VELGENDONE) + " and BOMDDONE = " + str(BOMDDONE))
        logging.info("PID" + str(os.getpid()) + ": Calc were load from " + expectedObjPath)
        if VELGENWRONG:
            Calc.removeSubCalc(filepath)
        else:
            Calc.updateSubCalc(filepath)
            if Calc.end <= 2:
                if args.config != None:
                    Calc.processConfig(args.config)
                    if Calc.end <= 2:
                        logging.warn("PID" + str(os.getpid()) + ": Session was originally set ending at point " + str(
                            Calc.end) + ". In order to continue, please specify a config file with end>2.")
                        logging.info("PID" + str(os.getpid()) + ": Exiting GSTA\n")
                        sys.exit(0)
                else:
                    logging.warn("PID" + str(os.getpid()) + ": Session was originally set ending at point " + str(
                        Calc.end) + ". In order to continue, please specify a config file with end>2.")
                    logging.info("PID" + str(os.getpid()) + ": Exiting GSTA\n")
                    sys.exit(0)
            if Calc.end > 2:
                if "_velgen." in filepath:
                    name = "traj" + filepath.split("_traj")[1].split('_velgen.')[0]
                else:
                    name = "traj" + filepath.split("_traj")[1].split('.')[0]
                for SC in Calc.SubCalcs:
                    if SC.name == name:
                        g.createMDInput(SC)
                        # if Calc.end > 3:
                        #     g.runMD_old(SC)
        Calc.save(lock=FL)
        if Calc.end > 3 and Calc.isVELGENDONE():
            g.runMD(Calc)
        logging.info("PID" + str(os.getpid()) + ": Exiting GSTA\n")
        sys.exit(0)

    elif MOVE:
        fullpath = os.path.realpath(args.file)
        tree = fullpath.split(os.sep)
        expectedObjPath = ""
        for i in range(len(tree) - 2): expectedObjPath += tree[i] + os.sep
        chks = glob.glob(expectedObjPath + os.sep + ".CalcFile_*")
        if len(chks) == 0:
            logging.warn("Pickled Calc File not found at " + expectedObjPath)
            logging.warn("Exiting")
            sys.exit(0)
        elif len(chks) == 1:
            expectedObjPath = chks[0]
        else:
            logging.warn(
                "More then one pickled Calc File found at " + expectedObjPath + ". Please remove those not linked to " + args.file)
            logging.warn("Exiting")
            sys.exit(0)
        try:
            Calc, FL = Calc.load(expectedObjPath)
        except IOError:
            sys.exit(0)

        os.chdir(Calc.workdir)
        logging.basicConfig(filename=Calc.freqfile.split('.')[0] + ".lg", level=logging.INFO,
                            format='%(asctime)s %(message)s')
        logging.info("PID" + str(os.getpid()) + ": " + programName + " started with MOVE = " + str(MOVE))
        logging.info("PID" + str(os.getpid()) + ": Calc were load from " + expectedObjPath)
        if VELGENWRONG:
            Calc.removeSubCalc(filepath)
        else:
            Calc.updateSubCalc(filepath)
        name = "traj" + filepath.split("_traj")[1].split('.')[0]
        for SC in Calc.SubCalcs:
            if SC.name == name:
                g.createMovedVelInputs(SC)
        Calc.save(lock=FL)

        if Calc.isMOVED():
            g.runVelGen(Calc)

        logging.info("PID" + str(os.getpid()) + ": Exiting GSTA\n")
        sys.exit(0)

    else:
        try:
            Calc.mol.vibfreqs
        except AttributeError:
            print("No vibrations found. Is it a freq output?")
            print("Check config file if you wish to continue a session.")
            sys.exit(0)
        os.chdir(Calc.workdir)
        logging.basicConfig(filename=Calc.freqfile.split('.')[0] + ".lg", level=logging.INFO,
                            format='%(asctime)s %(message)s')
        logging.info("PID" + str(os.getpid()) + ": " + programName + " started.")
        if args.config != None:
            Calc.processConfig(args.config)
        if Calc.command == None: Calc.setCommand()
        if Calc.temp == None: Calc.setTemp()
        # if Calc.type == None: Calc.setType()
        if Calc.end == None: Calc.setEnd()
        if Calc.NTraj == None: Calc.setNTraj()
        # if Calc.type == 1:
        #     if Calc.maxvib == None: Calc.setMaxVib()
        # elif Calc.type == 2:
        #     if Calc.NTraj == None: Calc.setNTraj()
        # if Calc.type == 1:
        #     if Calc.end > 0:
        #         g.createVelInputs(Calc)
        #     if Calc.end > 1:
        #         g.runVelGen(Calc)
        # elif Calc.type == 2:
        if not Calc.rot:
            if Calc.end > 0:
                g.createMoveInputs(Calc)
            if Calc.end > 1:
                Calc.save()
                g.runMove(Calc)
        else:
            if Calc.end > 0:
                g.createVelInputs(Calc)
            if Calc.end > 1:
                Calc.save()
                g.runVelGen(Calc)

    return


if __name__ == "__main__": main()
