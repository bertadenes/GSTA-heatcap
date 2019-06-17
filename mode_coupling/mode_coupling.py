import os
import numpy as np
import modules.g09BOMD_filter as g
import modules.GetHeatcapacity_util as u
from modules.Heatcapacity_subproc import DOS_from_velocity_single
import mode_coupling.vibspectrum as spect
import matplotlib
matplotlib.use('TkAgg') # MUST BE CALLED BEFORE IMPORTING plt
import matplotlib.pyplot as plt


class Namespace:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


def coupling_process(freqname):
    freq = g.processg09output(freqname)
    zero_pot = freq.scfenergies[-1] * 0.036749324
    # MDname = os.path.join(os.path.abspath(os.path.dirname(freqname)), "MD",
    #                       "{0:s}_traj{1:d}_MD_pos.out".format(freqname.split(os.sep)[-1].split('.')[0], 2))
    MDname = "/media/berta/fourier/BPDT/BOMD/bpdt_2Au_mode48kin3000K.out"
    dos = DOS_from_velocity_single(MDname)
    params = spect.Params()
    params.files = [freqname]
    raman = spect.process(params)
    IR = spect.process_IR(params)
    x1 = range(6000)
    f2, ax2 = plt.subplots()
    plt.xlabel("Frequency cm^-1")
    plt.ylabel("Normalized Intensity")
    ax2.plot(x1, dos/np.max(dos), label="VDoS")
    ax2.plot(x1, raman[0]/np.max(raman[0]), label="Raman")
    ax2.plot(x1, IR[0]/np.max(IR[0]), label="IR")
    plt.legend(loc='best')
    plt.show()

    traj = Namespace()
    traj.coord = g.getCoords(MDname)[0, :, :, :]
    traj.vel = g.getVel(MDname)
    traj.Epot = g.getEpot(MDname)[0, :]
    traj.Ekin = np.mean(np.array(u.calcEkin(traj.vel)), axis=(0, 1))[:, 0]
    traj.vel = traj.vel[0, :, :, :]
    modecoeff = np.empty(shape=(len(traj.Ekin), len(freq.vibfreqs)), dtype=np.float_)
    kincoeff = np.empty(shape=(len(traj.Ekin), len(freq.vibfreqs)), dtype=np.float_)
    potcoeff = np.empty(shape=(len(traj.Ekin), len(freq.vibfreqs)), dtype=np.float_)
    for i in range(len(traj.Ekin)):
        # print("{0:f} {1:f} {2:f}".format(traj.Ekin[i],
        #                                  traj.Epot[i] - zero_pot,
        #                                  traj.Ekin[i]+traj.Epot[i]))
        # disp = np.reshape((traj.coord[i] - freq.atomcoords[-1]), freq.natom * 3)
        disp = np.reshape((traj.coord[i] - traj.coord[0]), freq.natom * 3)
        proj_pot = []
        proj_kin = []
        for v in freq.vibdisps:
            proj_pot.append(np.abs(np.dot(np.reshape(v, freq.natom * 3) / np.linalg.norm(np.reshape(v, freq.natom * 3)), disp)))
            proj_kin.append(np.abs(np.dot(np.reshape(v, freq.natom * 3) / np.linalg.norm(np.reshape(v, freq.natom * 3)), np.reshape(traj.vel[i], freq.natom * 3))))
        kinpart = traj.Ekin[i] / (traj.Ekin[i] + traj.Epot[i] - zero_pot)
        potpart = 1 - kinpart
        for j in range(len(freq.vibdisps)):
            kincoeff[i][j] = proj_kin[j] / np.sum(proj_kin)
            potcoeff[i][j] = proj_pot[j] / np.sum(proj_pot)
            modecoeff[i][j] = kinpart * kincoeff[i][j] + potpart * potcoeff[i][j]
        print(i+1, "done")
    f, ax = plt.subplots()
    plt.xlabel("Step")
    plt.ylabel("Coefficient")
    ax.plot(modecoeff[:, 57], label="total")
    ax.plot(kincoeff[:, 57], label="kinetic")
    ax.plot(potcoeff[:, 57], label="potential")
    # for k in range(66):
        # ax.plot(modecoeff[:, k], label="mode {:d}".format(k + 1))
        # ax.plot(kincoeff[:, k], label="mode {:d}".format(k + 1))
        # ax.plot(potcoeff[1:, k], label="mode {:d}".format(k + 1))
    ind = 0
    for p in np.max(potcoeff[1:, :], axis=0):
        ind += 1
        print(ind, p)
    plt.legend(loc='best')
    plt.show()

    print("stop")
