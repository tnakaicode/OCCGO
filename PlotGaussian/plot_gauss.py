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

        nx, ny = [int(v) for v in getline(self.cfg_txt, 1).split()]
        xs, xe = [float(v) for v in getline(self.cfg_txt, 2).split()]
        ys, ye = [float(v) for v in getline(self.cfg_txt, 3).split()]
        px = np.linspace(-xs, xe, nx)
        py = np.linspace(-ys, ye, ny)
        self.mesh = np.meshgrid(px, py)

        self.setup_freq()
        self.ampl = self.setup_ampl()
        self.phas = self.setup_phas()
        self.func = self.ampl * np.exp(1j * self.knum * self.phas)

    def setup_freq(self):
        w_val, w_unit = getline(cfg_txt, 10).split()
        if w_unit in ["GHz", "kHz", "Hz"]:
            self.freq = float(w_val)
            self.knum, self.wave = knum_from_freq(self.freq, w_unit)
        elif w_unit in ["m", "cm", "mm", "um"]:
            self.wave = float(w_val)
            self.knum, self.freq = knum_from_wave(self.wave, w_unit)
        else:
            self.freq = 100.0
            self.knum, self.wave = knum_from_freq(self.freq, "GHz")

    def setup_ampl(self):
        val = [float(v) for v in getline(self.cfg_txt, 5).split()]
        wxy = [float(v) for v in getline(self.cfg_txt, 6).split()]
        sxy = [float(v) for v in getline(self.cfg_txt, 7).split()]
        deg = float(getline(self.cfg_txt, 8).split()[0])
        ampl = gauss_2d(self.mesh, sxy, wxy, np.rad2deg(deg)) * val[0] + val[1]
        return ampl

    def setup_phas(self):
        p_prof = getline(self.cfg_txt, 11).split()[0]
        p_txy = [float(v) for v in getline(self.cfg_txt, 12).split()]
        p_sxy = [float(v) for v in getline(self.cfg_txt, 13).split()]
        if p_prof == "curvature":
            phas_x = curvature(self.mesh[0], p_txy[0], p_sxy[0])
            phas_y = curvature(self.mesh[1], p_txy[1], p_sxy[1])
        elif p_prof == "tilt":
            phas_x = self.mesh[0] * np.tan(np.deg2rad(p_txy[0]))
            phas_y = self.mesh[1] * np.tan(np.deg2rad(p_txy[1]))
        else:
            phas_x = np.zeros_like(self.mesh[0])
            phas_y = np.zeros_like(self.mesh[0])
        phas = phas_x + phas_y
        return phas


if __name__ == '__main__':
    argvs = sys.argv
    parser = OptionParser()
    parser.add_option("--dir", dest="dir", default="./")
    parser.add_option("--file", dest="file", default="plot_gauss.txt")
    opt, argc = parser.parse_args(argvs)
    print(opt, argc)

    cfg_txt = opt.file
    obj = GaussianProfile(cfg_txt)
    obj.contourf_div(obj.mesh, np.abs(obj.func), pngname=obj.tempname)
