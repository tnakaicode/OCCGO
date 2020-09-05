import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
from linecache import getline, clearcache
from optparse import OptionParser

sys.path.append(os.path.join("../"))
from src.base import plot2d, plot3d

import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)


def gauss_1d(px, sx=0, wx=10):
    py = np.exp(-0.5 * ((px - sx) / wx)**2)
    return py


def gauss_2d(mesh, sxy=[0, 0], wxy=[10, 10], deg=0.0):
    x, y = mesh[0] - sxy[0], mesh[1] - sxy[1]
    rot = np.deg2rad(deg)
    px = x * np.cos(rot) - y * np.sin(rot)
    py = y * np.cos(rot) + x * np.sin(rot)
    fx = np.exp(-0.5 * (px / wxy[0])**2)
    fy = np.exp(-0.5 * (py / wxy[1])**2)
    return fx * fy


if __name__ == '__main__':
    argvs = sys.argv
    parser = OptionParser()
    parser.add_option("--dir", dest="dir", default="./")
    parser.add_option("--pxyz", dest="pxyz",
                      default=[0.0, 0.0, 0.0], type="float", nargs=3)
    opt, argc = parser.parse_args(argvs)
    print(opt, argc)

    cfg_txt = "plot_gauss.txt"
    lx, ly = [float(v) for v in getline(cfg_txt, 1).split()]
    sx, sy = [float(v) for v in getline(cfg_txt, 2).split()]
    nx, ny = [int(v) for v in getline(cfg_txt, 3).split()]
    px = np.linspace(-1, 1, nx) * lx / 2 - sx
    py = np.linspace(-1, 1, ny) * ly / 2 - sy
    mesh = np.meshgrid(px, py)

    val, = [float(v) for v in getline(cfg_txt, 5).split()]
    wxy = [float(v) for v in getline(cfg_txt, 6).split()]
    sxy = [float(v) for v in getline(cfg_txt, 7).split()]
    deg, = [float(v) for v in getline(cfg_txt, 8).split()]
    func = gauss_2d(mesh, sxy, wxy, np.rad2deg(deg)) * val

    obj = plot2d()
    obj.contourf_div(mesh, func, pngname=obj.tempname)
