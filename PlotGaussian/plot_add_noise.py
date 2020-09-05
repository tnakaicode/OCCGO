import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
from linecache import getline, clearcache
from optparse import OptionParser

from sympy import im

sys.path.append(os.path.join("../"))
from src.base import plot2d, plot3d
from src.profile import gauss_1d, gauss_1d_skew

import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)


def make_noise(px):
    py = np.zeros_like(px)
    for i, x in enumerate(px):
        py[i] = np.random.random() * (1 - x)**2 / 5
    return py


if __name__ == '__main__':
    argvs = sys.argv
    parser = OptionParser()
    parser.add_option("--dir", dest="dir", default="./")
    parser.add_option("--pxyz", dest="pxyz",
                      default=[0.0, 0.0, 0.0], type="float", nargs=3)
    opt, argc = parser.parse_args(argvs)
    print(opt, argc)

    px = np.linspace(-1, 1, 100)
    py = gauss_1d(px, wx=0.1) + 0.1
    #py = py - py.min()
    py = py / py.max()
    obj = plot2d(aspect="equal")
    obj.axs.set_ylim(0,1.1)
    obj.axs.plot(px, py)
    obj.axs.plot(px, py + make_noise(py))
    obj.SavePng()
