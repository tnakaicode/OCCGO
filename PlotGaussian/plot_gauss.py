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
from src.Unit import knum_from_freq, knum_from_wave
from src.profile import gaussian_func, rot_mesh
from src.geomtory import curvature
from src.Gaussian import GaussianProfile

import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)


def gen_noise(px=[0, 0], ratio=100):
    print(px.shape)
    py = np.random.random(px.shape)
    py = py * 2 - 1
    py = py * (1 - px)**2 * ratio
    return px + py


if __name__ == '__main__':
    argvs = sys.argv
    parser = OptionParser()
    parser.add_option("--dir", dest="dir", default="./")
    parser.add_option("--file", dest="file", default="plot_gauss.txt")
    parser.add_option("--rati", dest="rati", default=1.0E-03, type="float")
    opt, argc = parser.parse_args(argvs)
    print(opt, argc)

    cfg_txt = opt.file
    ratio = opt.rati
    obj = GaussianProfile(cfg_txt)
    obj.ampl = gen_noise(obj.ampl, ratio)
    obj.func = obj.ampl * np.exp(1j * obj.phas)
    obj.g_func = obj.create_gauss()
    obj.tmpdir = "./img/"
    obj.tempname = obj.tmpdir + "gaussian"
    obj.profile_out()
