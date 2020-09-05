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


if __name__ == '__main__':
    argvs = sys.argv
    parser = OptionParser()
    parser.add_option("--dir", dest="dir", default="./")
    parser.add_option("--file", dest="file", default="plot_gauss.txt")
    opt, argc = parser.parse_args(argvs)
    print(opt, argc)

    cfg_txt = opt.file
    obj = GaussianProfile(cfg_txt)
    obj.tmpdir = "./img/"
    obj.tempname = obj.tmpdir + "gaussian"
    obj.profile_out()
