from OCC.Core.gp import gp_Parab
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
from src.profile import gauss_2d, gauss_2d_skew

import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)


def gen_noise(px=[0, 0], ratio=100):
    print(px.shape)
    py = np.random.random(px.shape)
    py = py * 2 - 1
    py = py * (1 - px)**2 / ratio
    return px + py


if __name__ == '__main__':
    argvs = sys.argv
    parser = OptionParser()
    parser.add_option("--dir", dest="dir", default="./")
    parser.add_option("--pxyz", dest="pxyz",
                      default=[0.0, 0.0, 0.0], type="float", nargs=3)
    opt, argc = parser.parse_args(argvs)
    print(opt, argc)

    px = np.linspace(-1, 1, 100)
    mesh = np.meshgrid(px, px)
    pg_base = 5.0 * 10**(-4)
    pg_1d = gauss_1d(px, wx=0.1) + 0.1
    pg_1d = pg_1d / pg_1d.max()
    pg_2d = gauss_2d(mesh, wxy=[0.1, 0.1]) + pg_base
    pg_nd = gen_noise(pg_2d, 0.9 / pg_base)
    pg_db = 10 * np.log10(pg_nd)

    obj = plot2d(aspect="equal")
    obj.axs.set_ylim(0, 1.1)
    obj.axs.plot(px, pg_1d)
    obj.axs.plot(px, gen_noise(pg_1d, 1.0 / 0.1))
    obj.SavePng(obj.tempname + "_1d.png")

    obj.contourf_div(mesh, pg_2d, pngname=obj.tempname + "_2d.png")
    obj.contourf_div(mesh, pg_nd, pngname=obj.tempname + "_2d_noise.png")
    obj.contourf_div(mesh, pg_db, pngname=obj.tempname + "_2d_10db.png")
