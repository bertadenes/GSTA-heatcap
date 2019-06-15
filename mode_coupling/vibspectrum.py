from tkinter import *
from tkinter import filedialog
import cclib
import numpy as np
from scipy import constants
import matplotlib
matplotlib.use('TkAgg') # MUST BE CALLED BEFORE IMPORTING plt
import matplotlib.pyplot as plt
import os


class Params:
    def __init__(self):
        self.temperature = 298.15
        self.min = 0
        self.max = 5999
        self.sigma = 5.0
        self.source = 633.0
        self.files = []

    def update(self, temp, min, max, sigma, source):
        self.temperature = float(temp)
        self.min = int(min)
        self.max = int(max)
        self.sigma = float(sigma)
        self.source = float(source)

    def openFiles(self):
        # Create Tk root
        root = Tk()
        # Hide the main window
        root.withdraw()
        root.call('wm', 'attributes', '.', '-topmost', True)
        self.files = filedialog.askopenfilename(multiple=True)


def gauss(wn, x, sigma):
    y = np.exp(-2*(np.square(wn-x)/np.square(sigma)))
    return y


def getBroad(spectralmax, sigma, freq, raman):
    curve = np.zeros(shape=(spectralmax+1), dtype=np.float_)
    wn = np.arange((spectralmax+1))
    for i in range(len(freq)):
        curve = curve + raman[i]*gauss(wn, freq[i], sigma=sigma)
    return curve


def process(args):
    broads = []
    for file in args.files:
        data = cclib.io.ccread(file)
        intensity = (np.power(2 * constants.pi * ((1e7 / args.source) - data.vibfreqs), 4) / 45) * (
                constants.h / (8 * (constants.pi ** 2) * 100 * constants.c * data.vibfreqs * (
                1 - np.exp(- (constants.h * data.vibfreqs * 100 * constants.c) / (constants.k * args.temperature))
        ))) * (data.vibramans / constants.atomic_mass)
        broads.append(getBroad(args.max, args.sigma, data.vibfreqs, intensity)[args.min:(args.max + 1)])
    return broads


def process_IR(args):
    broads = []
    for file in args.files:
        data = cclib.io.ccread(file)
        intensity = data.vibirs
        broads.append(getBroad(args.max, args.sigma, data.vibfreqs, intensity)[args.min:(args.max + 1)])
    return broads


def plot(args, broads):
    wn = np.arange(args.min, args.max + 1)
    f, ax = plt.subplots()
    plt.xlabel("Wavenumber cm^-1")
    plt.ylabel("Raman intensity")
    plt.axis([args.min, args.max + 1, 0, np.max(broads)])
    for i in range(len(broads)):
        ax.plot(wn, broads[i], label=args.files[i].split(os.sep)[-1])
    plt.legend(loc='upper left')
    plt.show()


def main():
    return


if __name__ == '__main__':
    main()
