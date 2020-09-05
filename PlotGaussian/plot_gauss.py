import shutil
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
    p_x = px / px.max()
    py = np.random.random(px.shape)
    py = py * 2 - 1
    py = py * (1 - p_x)**2 * ratio
    return px + py * px.max()


def copy_file(str1, str2):
    if os.path.exists(str2):
        pass
    else:
        shutil.copyfile(str1, str2)


if __name__ == '__main__':
    argvs = sys.argv
    parser = OptionParser()
    parser.add_option("--dir", dest="dir", default="./img/")
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

    obj.tmpdir = "."
    obj.tmpdir = obj.add_dir(opt.dir)
    obj.tempname = obj.tmpdir + "gaussian"
    obj.profile_out()

    fp = open(obj.tmpdir + "noise.txt", "w")
    fp.write("{:.5E}\n".format(ratio))
    fp.close()
    copy_file(obj.cfg_txt, obj.tempname + ".txt")
    copy_file("./img/NOTE.md", obj.tmpdir + "NOTE.md")
