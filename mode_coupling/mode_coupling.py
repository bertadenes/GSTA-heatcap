import os
import numpy as np
import modules.g09BOMD_filter as g
import modules.GetHeatcapacity_util as u
import matplotlib
matplotlib.use('TkAgg') # MUST BE CALLED BEFORE IMPORTING plt
import matplotlib.pyplot as plt


class Namespace:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


def coupling_process(freqname):
    freq = g.processg09output(freqname)
    zero_pot = freq.scfenergies[-1] * 0.036749324
    MDname = os.path.join(os.path.abspath(os.path.dirname(freqname)), "MD",
                          "{0:s}_traj{1:d}_MD_pos.out".format(freqname.split(os.sep)[-1].split('.')[0], 1))
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
        disp = np.reshape((traj.coord[i] - freq.atomcoords[-1]), freq.natom * 3)
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
    plt.legend(loc='upper left')
    plt.show()

    print("stop")
