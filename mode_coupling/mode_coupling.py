import os
import numpy as np
import modules.g09BOMD_filter as g
import modules.GetHeatcapacity_util as u


class Namespace:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


def coupling_process(freqname):
    freq = g.processg09output(freqname)
    MDname = os.path.join(os.path.abspath(os.path.dirname(freqname)), "MD",
                          "{0:s}_traj{1:d}_MD_pos.out".format(freqname.split(os.sep)[-1].split('.')[0], 1))
    traj = Namespace()
    traj.coord = g.getCoords(MDname)[0, :, :, :]
    traj.vel = g.getVel(MDname)
    traj.Epot = g.getEpot(MDname)[0, :]
    traj.Ekin = np.mean(np.array(u.calcEkin(traj.vel)), axis=(0, 1))[:, 0]
    traj.vel = traj.vel[0, :, :, :]
    for i in range(len(traj.Ekin)):
        print("{0:f} {1:f} {2:f}".format(traj.Ekin[i],
                                         traj.Epot[i], #freq.scfenergies[-1]) * 0.036749324,
                                         traj.Ekin[i]+traj.Epot[i]))

    print("stop")