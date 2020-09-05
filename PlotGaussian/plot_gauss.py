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
from src.geomtory import curvature

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


class GaussianProfile(plot2d):

    def __init__(self, cfgtxt="plot_gauss.txt"):
        plot2d.__init__(self)
        self.cfg_txt = cfgtxt

        lx, ly = [float(v) for v in getline(self.cfg_txt, 1).split()]
        sx, sy = [float(v) for v in getline(self.cfg_txt, 2).split()]
        nx, ny = [int(v) for v in getline(self.cfg_txt, 3).split()]
        px = np.linspace(-1, 1, nx) * lx / 2 - sx
        py = np.linspace(-1, 1, ny) * ly / 2 - sy
        self.mesh = np.meshgrid(px, py)
        self.ampl = self.setup_ampl()

    def setup_ampl(self):
        val = [float(v) for v in getline(self.cfg_txt, 5).split()]
        wxy = [float(v) for v in getline(self.cfg_txt, 6).split()]
        sxy = [float(v) for v in getline(self.cfg_txt, 7).split()]
        deg = float(getline(self.cfg_txt, 8).split()[0])
        ampl = gauss_2d(mesh, sxy, wxy, np.rad2deg(deg)) * val[0] + val[1]
        return ampl


if __name__ == '__main__':
    argvs = sys.argv
    parser = OptionParser()
    parser.add_option("--dir", dest="dir", default="./")
    parser.add_option("--file", dest="file", default="plot_gauss.txt")
    opt, argc = parser.parse_args(argvs)
    print(opt, argc)

    cfg_txt = opt.file
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
    ampl = gauss_2d(mesh, sxy, wxy, np.rad2deg(deg)) * val + 1.0

    w_val, w_unit = getline(cfg_txt, 10).split()
    if w_unit in ["GHz", "kHz", "Hz"]:
        freq = float(w_val)
        knum, wave = knum_from_freq(freq, w_unit)
    elif w_unit in ["m", "cm", "mm", "um"]:
        wave = float(w_val)
        knum, freq = knum_from_wave(wave, w_unit)
    else:
        freq = 100.0
        knum, wave = knum_from_freq(freq, "GHz")

    p_prof = getline(cfg_txt, 11).split()[0]
    p_txy = [float(v) for v in getline(cfg_txt, 12).split()]
    p_sxy = [float(v) for v in getline(cfg_txt, 13).split()]
    if p_prof == "curvature":
        phas_x = curvature(mesh[0], p_txy[0], p_sxy[0])
        phas_y = curvature(mesh[1], p_txy[1], p_sxy[1])
    elif p_prof == "tilt":
        phas_x = mesh[0] * np.tan(np.deg2rad(p_txy[0]))
        phas_y = mesh[1] * np.tan(np.deg2rad(p_txy[1]))
    else:
        phas_x = np.zeros_like(ampl)
        phas_y = np.zeros_like(ampl)
    phas = phas_x + phas_y
    func = ampl * np.exp(1j * knum * phas)

    obj = plot2d()
    obj.contourf_div(mesh, np.abs(func), pngname=obj.tempname)
