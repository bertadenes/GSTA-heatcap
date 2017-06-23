#!/usr/bin/env python3
"""Generic linux daemon base class for python 3.x."""

import sys, os, time, atexit, signal, glob, logging, argparse
from subprocess import Popen

class Daemon:
    """A generic daemon class.
    Usage: subclass the daemon class and override the run() method."""

    def __init__(self, pidfile): self.pidfile = pidfile
    
    def daemonize(self):
        """Deamonize class. UNIX double fork mechanism."""

        try: 
            pid = os.fork() 
            if pid > 0:
                # exit first parent
                sys.exit(0) 
        except OSError as err: 
            sys.stderr.write('fork #1 failed: {0}\n'.format(err))
            sys.exit(1)
    
        # decouple from parent environment
        os.chdir('/') 
        os.setsid() 
        os.umask(0) 
    
        # do second fork
        try: 
            pid = os.fork() 
            if pid > 0:

                # exit from second parent
                sys.exit(0) 
        except OSError as err: 
            sys.stderr.write('fork #2 failed: {0}\n'.format(err))
            sys.exit(1) 
    
        # redirect standard file descriptors
        sys.stdout.flush()
        sys.stderr.flush()
        si = open(os.devnull, 'r')
        so = open(os.devnull, 'a+')
        se = open(os.devnull, 'a+')

        os.dup2(si.fileno(), sys.stdin.fileno())
        os.dup2(so.fileno(), sys.stdout.fileno())
        os.dup2(se.fileno(), sys.stderr.fileno())
    
        # write pidfile
        atexit.register(self.delpid)

        pid = str(os.getpid())
        with open(self.pidfile,'w+') as f:
            f.write(pid + '\n')
    
    def delpid(self):
        os.remove(self.pidfile)

    def start(self):
        """Start the daemon."""

        # Check for a pidfile to see if the daemon already runs
        try:
            with open(self.pidfile,'r') as pf:

                pid = int(pf.read().strip())
        except IOError:
            pid = None
    
        if pid:
            message = "pidfile {0} already exist. " + \
                    "Daemon already running?\n"
            sys.stderr.write(message.format(self.pidfile))
            sys.exit(1)
        
        # Start the daemon
        self.daemonize()
        self.run()

    def stop(self):
        """Stop the daemon."""

        # Get the pid from the pidfile
        try:
            with open(self.pidfile,'r') as pf:
                pid = int(pf.read().strip())
        except IOError:
            pid = None
    
        if not pid:
            message = "pidfile {0} does not exist. " + \
                    "Daemon not running?\n"
            sys.stderr.write(message.format(self.pidfile))
            return # not an error in a restart

        # Try killing the daemon process    
        try:
            while 1:
                os.kill(pid, signal.SIGTERM)
                time.sleep(0.1)
        except OSError as err:
            e = str(err.args)
            if e.find("No such process") > 0:
                if os.path.exists(self.pidfile):
                    os.remove(self.pidfile)
            else:
                print (str(err.args))
                sys.exit(1)

    def restart(self):
        """Restart the daemon."""
        self.stop()
        self.start()

    def run(self):
        """You should override this method when you subclass Daemon.
        
        It will be called after the process has been daemonized by 
        start() or restart()."""

def tail(f,BLOCK_SIZE = 256):
    f.seek(0,2)
    size = f.tell()
    if size < BLOCK_SIZE: BLOCK_SIZE = size
    f.seek(-BLOCK_SIZE, 2)
    a = f.readlines()
    rv = ""
    for b in a: rv += str(b)
    return rv

class g09daemon(Daemon):
    """created 25/07/2016
    deamon for watching g09 output"""

    def __init__(self,pidfile,dirname,inputname,code=None):
        self.pidfile = pidfile
        self.dirname = dirname
        self.inputname = inputname
        self.code = code

    def run(self):
        """doc"""
        self.filename = self.inputname.split('.')[0]
        self.minuteCount = 0
        self.outpath = None
        logging.basicConfig(filename=self.dirname+os.sep+"g09daemon.lg",level=logging.INFO, format='%(asctime)s %(message)s')
        logging.info(self.filename+": Daemon started. PID: "+str(os.getpid()))
        while True:
            if len(glob.glob(self.dirname+os.sep+self.filename+"*log")) == 1 :
                self.outpath = glob.glob(self.dirname+os.sep+self.filename+"*log")[0]
                logging.info(self.filename+": Output found at "+self.outpath)#should write into a log
                time.sleep(10)
                break
            elif len(glob.glob(self.dirname+os.sep+self.filename+"*out")) == 1:
                self.outpath = glob.glob(self.dirname+os.sep+self.filename+"*out")[0]
                logging.info(self.filename+": Output found at "+self.outpath)#should write into a log
                time.sleep(10)
                break
            elif self.minuteCount == 720:
                #output does not show up for 12 hours
                logging.info(self.filename+": Output have not been found for 12 hours. Exiting")
                sys.exit(0)
            else:
                time.sleep(60)
                self.minuteCount += 1
        while True:
            with open(self.outpath,'rb') as outfile:
                if "termination" in tail(outfile):
                    if "Error" in tail(outfile):
                        logging.info(self.filename+": Calculation terminated with error.")
                        sys.exit(0)
                    else:
                        if self.code == None:
                            logging.info(self.filename+": Normal termination. Recalling main program: "+"entropy -f "+str(self.outpath))
                            Popen(["entropy","-f",str(self.outpath)])
                        else:
                            logging.info(self.filename+": Normal termination. Recalling main program: "+"entropy -f "+str(self.outpath)+" -t "+str(self.code))
                            Popen(["entropy","-f",str(self.outpath),"-t",str(self.code)])
                        logging.info(self.filename+": Exiting. PID: "+str(os.getpid()))
                        sys.exit(0)
                else:
                    time.sleep(60)
        return

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('dir')
    parser.add_argument('inp')
    parser.add_argument('-t','--type')
    args = parser.parse_args()
    daemon = g09daemon(pidfile=None,dirname=args.dir,inputname=args.inp,code=args.type)
    g09daemon.run(daemon)

if __name__ == "__main__": main()
