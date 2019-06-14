#!/usr/bin/env python3

import glob
import logging
import os
import sys
from subprocess import Popen
from time import sleep

from daemon.g09daemon import tail


class NoOutputException(Exception):
    pass

class StillRunningException(Exception):
    pass

program = "GSTA-hc"

class g09calcHandler():
    def submit(self,inp):
        os.chdir(self.workdir)
        Popen([str(self.command), inp, self.options])
        #logging.info("PID" + str(os.getpid()) + ": " + inp + " was submitted")
        self.submittedInputs.append(inp)
        self.inputs.remove(inp)

    def seekOut(self,sinp):
        filename = sinp.split('.')[0]
        if len(glob.glob(os.path.join(self.workdir,filename) + ".log")) == 1:
            self.foundOutputs.append(glob.glob(os.path.join(self.workdir,filename) + "*log")[0])
            self.toBeRemoved.append(sinp)
            return
        elif len(glob.glob(os.path.join(self.workdir, filename) + ".out")) == 1:
            self.foundOutputs.append(glob.glob(os.path.join(self.workdir, filename) + "*out")[0])
            self.toBeRemoved.append(sinp)
            return
        else:
            raise NoOutputException

    def checkOut(self,fout):
        with open(fout, 'rb') as outfile:
            if "termination" in tail(outfile):
                if "Normal" in tail(outfile):
                    self.termOutputs.append(fout)
                    self.toBeRemoved.append(fout)
                    return
                if "Error" in tail(outfile):
                    self.wrongOutputs.append(fout)
                    return
                else:
                    logging.info("PID" + str(os.getpid()) + ": "+fout+" stopped unexpectedly.")
                    self.wrongOutputs.append(fout)
                    return
            else:
                raise StillRunningException

    def recallMain(self,tout):
        Popen([program, "-f", tout, "-t", "velgen"])
        self.recalledOutputs.append(tout)
        self.toBeRemoved.append(tout)

    def watch(self):
        for inp in self.inputs:
            try:
                self.seekOut(inp)
            except NoOutputException:
                pass
        for rem in self.toBeRemoved:
            self.inputs.remove(rem)
        self.toBeRemoved = []
        for fout in self.foundOutputs:
            try:
                self.checkOut(fout)
            except StillRunningException:
                pass
        for rem in self.toBeRemoved:
            self.foundOutputs.remove(rem)
        self.toBeRemoved = []
        if self.RECALL:
            for tout in self.termOutputs:
                self.recallMain(tout)
                sleep(10)
            for rem in self.toBeRemoved:
                self.termOutputs.remove(rem)
            self.toBeRemoved = []
        logging.info("PID" + str(os.getpid()) + ": Watch started. The current state of calculations:")
        logging.info("PID" + str(os.getpid()) + ": total calculations: " + str(self.all))
        logging.info("PID" + str(os.getpid()) + ": submitted: " + str(len(self.submittedInputs)))
        logging.info("PID" + str(os.getpid()) + ": running: " + str(len(self.foundOutputs)))
        logging.info("PID" + str(os.getpid()) + ": normally terminated: " + str(len(self.termOutputs)))
        logging.info("PID" + str(os.getpid()) + ": sent back: " + str(len(self.recalledOutputs)))
        logging.info("PID" + str(os.getpid()) + ": terminated with error: " + str(len(self.wrongOutputs)))
        logging.info("PID" + str(os.getpid()) + ": remaining: " + str(len(self.inputs)))
        if self.RECALL:
            done = len(self.wrongOutputs) + len(self.recalledOutputs)
        else:
            done = len(self.wrongOutputs) + len(self.termOutputs)
        cycles = 0
        while done != self.all:
            running = len(self.submittedInputs) + len(self.foundOutputs)
            while running < self.NCalc:
                if len(self.inputs) == 0:
                    break
                else:
                    self.submit(self.inputs[0])
                    running = len(self.submittedInputs) + len(self.foundOutputs)
            sleep(15)
            for sinp in self.submittedInputs:
                try:
                    self.seekOut(sinp)
                except NoOutputException:
                    pass
            for rem in self.toBeRemoved:
                self.submittedInputs.remove(rem)
            self.toBeRemoved = []
            sleep(15)
            for fout in self.foundOutputs:
                try:
                    self.checkOut(fout)
                except StillRunningException:
                    pass
            for rem in self.toBeRemoved:
                self.foundOutputs.remove(rem)
            self.toBeRemoved = []
            if self.RECALL:
                for tout in self.termOutputs:
                    self.recallMain(tout)
                    sleep(10)
                for rem in self.toBeRemoved:
                    self.termOutputs.remove(rem)
                self.toBeRemoved = []
            cycles += 1
            if cycles == 40:
                logging.info("PID" + str(os.getpid()) + ": total calculations: " + str(self.all))
                logging.info("PID" + str(os.getpid()) + ": submitted: " + str(len(self.submittedInputs)))
                logging.info("PID" + str(os.getpid()) + ": running: " + str(len(self.foundOutputs)))
                logging.info("PID" + str(os.getpid()) + ": normally terminated: " + str(len(self.termOutputs)))
                logging.info("PID" + str(os.getpid()) + ": sent back: " + str(len(self.recalledOutputs)))
                logging.info("PID" + str(os.getpid()) + ": terminated with error: " + str(len(self.wrongOutputs)))
                logging.info("PID" + str(os.getpid()) + ": remaining: " + str(len(self.inputs)))
                cycles = 0
            if self.RECALL:
                done = len(self.wrongOutputs) + len(self.recalledOutputs)
            else:
                done = len(self.wrongOutputs) + len(self.termOutputs)
        logging.info("\nPID" + str(os.getpid()) + ": The shift is over, watchman may rest.")
        logging.info("PID" + str(os.getpid()) + ": total calculations: " + str(self.all))
        logging.info("PID" + str(os.getpid()) + ": sent back: " + str(len(self.recalledOutputs)))
        logging.info("PID" + str(os.getpid()) + ": not sent back: " + str(len(self.termOutputs)))
        logging.info("PID" + str(os.getpid()) + ": terminated with error: " + str(len(self.wrongOutputs))+"\n")
        if len(self.wrongOutputs) != 0:
            for wout in self.wrongOutputs:
                Popen([program, "-f", wout, "-t", "velgenwrong"])
        sys.exit(0)

class g09velgenHandler(g09calcHandler):
    def __init__(self,Calc):
        if len(Calc.velGenInputs) != 0:
            self.inputs = Calc.velGenInputs
        else:
            self.inputs = []
            for SC in Calc.SubCalcs:
                self.inputs.append(SC.velgenInput)
        self.all = len(self.inputs)
        self.submittedInputs = []
        self.foundOutputs = []
        self.termOutputs = []
        self.recalledOutputs = []
        self.wrongOutputs = []
        self.toBeRemoved = []
        self.workdir = os.path.join(Calc.workdir, "velgen")
        self.command = Calc.command
        self.options = Calc.options
        self.NCalc = Calc.NCalc
        if Calc.end > 2:
            self.RECALL = True
        else: self.RECALL = False
        logging.basicConfig(filename=os.path.join(self.workdir,"g09daemon.lg"), level=logging.INFO,
                            format='%(asctime)s %(message)s')

class g09MDHandler(g09calcHandler):
    def __init__(self,Calc):
        self.inputs = []
        for SC in Calc.SubCalcs:
            self.inputs.append(SC.MDinput1)
            self.inputs.append(SC.MDinput2)
        self.all = len(self.inputs)
        self.submittedInputs = []
        self.foundOutputs = []
        self.termOutputs = []
        self.recalledOutputs = []
        self.wrongOutputs = []
        self.toBeRemoved = []
        self.workdir = os.path.join(Calc.workdir, "MD")
        self.command = Calc.command
        self.options = Calc.options
        self.NCalc = Calc.NCalc
        if Calc.end > 4:
            self.RECALL = True
        else: self.RECALL = False
        logging.basicConfig(filename=os.path.join(self.workdir, "g09daemon.lg"), level=logging.INFO,
                            format='%(asctime)s %(message)s')
    def recallMain(self,tout):
        Popen([program, "-f", tout, "-t", "BOMD"])
        self.recalledOutputs.append(tout)
        self.toBeRemoved.append(tout)

class g09MoveHandler(g09calcHandler):
    def __init__(self,Calc):
        self.inputs = Calc.moveInputs
        self.all = len(self.inputs)
        self.submittedInputs = []
        self.foundOutputs = []
        self.termOutputs = []
        self.recalledOutputs = []
        self.wrongOutputs = []
        self.toBeRemoved = []
        self.workdir = os.path.join(Calc.workdir, "velgen")
        self.command = Calc.command
        self.options = Calc.options
        self.NCalc = Calc.NCalc
        if Calc.end > 1:
            self.RECALL = True
        else: self.RECALL = False
        logging.basicConfig(filename=os.path.join(self.workdir, "g09daemon.lg"), level=logging.INFO,
                            format='%(asctime)s %(message)s')
    def recallMain(self,tout):
        Popen([program, "-f", tout, "-t", "move"])
        self.recalledOutputs.append(tout)
        self.toBeRemoved.append(tout)