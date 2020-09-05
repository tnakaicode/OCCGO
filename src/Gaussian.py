import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
from linecache import getline, clearcache
from optparse import OptionParser
from numpy.core.defchararray import array
from numpy.lib.function_base import unwrap

from sympy import im

sys.path.append(os.path.join("../"))
from src.base import plot2d, plot3d
from src.Unit import knum_from_freq, knum_from_wave

from src.geomtory import curvature
from src.profile import gaussian_func, rot_mesh
from src.profile import get_centroid, get_wxy

import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)


class GaussianProfile(plot2d):

    def __init__(self, cfgtxt="plot_gauss.txt"):
        plot2d.__init__(self)
        self.cfg_txt = cfgtxt
        self.setup_mesh()
        self.setup_freq()

        self.ampl = self.setup_ampl()
        self.phas = self.setup_phas()
        self.func = self.ampl * np.exp(1j * self.knum * self.phas)
        self.g_func = self.create_gauss()

    def setup_mesh(self):
        nx, ny = [int(v) for v in getline(self.cfg_txt, 1).split()]
        xs, xe = [float(v) for v in getline(self.cfg_txt, 2).split()]
        ys, ye = [float(v) for v in getline(self.cfg_txt, 3).split()]
        px = np.linspace(xs, xe, nx)
        py = np.linspace(ys, ye, ny)
        self.mesh = np.meshgrid(px, py)

    def setup_func(self):
        self.setup_freq()
        self.ampl = self.setup_ampl()
        self.phas = self.setup_phas()
        self.func = self.ampl * np.exp(1j * self.knum * self.phas)

    def setup_freq(self):
        w_val, w_unit = getline(self.cfg_txt, 10).split()
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
        ampl = gaussian_func(self.mesh, sxy, wxy, np.deg2rad(deg))
        ampl = ampl * val[0] + val[1]
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

    def create_gauss(self):
        ampl = np.abs(self.func)
        sxy = get_centroid(self.mesh, ampl)
        wxy, rot = get_wxy(self.mesh, ampl)
        g_ampl = gaussian_func(self.mesh, sxy, wxy, rot)
        g_func = g_ampl * np.exp(1j * np.zeros_like(g_ampl))
        return g_func

    def contourf_comp(self, mesh, func1, func2, loc=[0, 0], txt="", title="name", pngname=None, level=None):
        sx, sy = loc
        nx, ny = func1.shape
        xs, ys = mesh[0][0, 0], mesh[1][0, 0]
        xe, ye = mesh[0][0, -1], mesh[1][-1, 0]
        dx, dy = mesh[0][0, 1] - mesh[0][0, 0], mesh[1][1, 0] - mesh[1][0, 0]
        mx, my = int((sy - ys) / dy), int((sx - xs) / dx)
        tx, ty = 1.1, 0.0

        self.new_2Dfig()
        self.div_axs()
        self.ax_x.plot(mesh[0][mx, :], func1[mx, :])
        self.ax_x.plot(mesh[0][mx, :], func2[mx, :])
        self.ax_x.set_title("y = {:.2f}".format(sy))

        self.ax_y.plot(func1[:, my], mesh[1][:, my])
        self.ax_y.plot(func2[:, my], mesh[1][:, my])
        self.ax_y.set_title("x = {:.2f}".format(sx))

        if type(level) == np.ndarray:
            dl = level[2] - level[0]
            ls = level[0] - dl * 1.25
            le = level[-1] + dl * 1.25
            self.ax_x.set_ylim(ls, le)
            self.ax_y.set_xlim(ls, le)

        self.fig.text(tx, ty, txt, transform=self.ax_x.transAxes)
        im = self.axs.contourf(*mesh, func1, cmap="jet", levels=level)
        self.axs.set_title(title)
        self.fig.colorbar(im, ax=self.axs, shrink=0.9)

        plt.tight_layout()
        if pngname == None:
            self.SavePng_Serial(pngname)
        else:
            self.SavePng(pngname)

    def profile_out(self):
        ampl = np.abs(self.func)**2
        ampl_norm = ampl / ampl.max()
        db10 = 10 * np.log10(ampl_norm)
        phas = np.angle(self.func)
        phas_norm = unwrap(phas)
        sxy = get_centroid(self.mesh, ampl)
        wxy, rot = get_wxy(self.mesh, ampl)

        g_ampl = np.abs(self.g_func)**2
        g_ampl = g_ampl / g_ampl.max()
        g_db10 = 10 * np.log10(g_ampl)
        level = np.linspace(-50, 0, 11)

        name = self.tempname
        self.contourf_div(self.mesh, ampl, sxy,
                          pngname=name + "_ampl.png")
        self.contourf_div(self.mesh, db10, sxy,
                          pngname=name + "_10db.png", level=level)
        self.contourf_div(self.mesh, phas, sxy,
                          pngname=name + "_phas.png")
        self.contourf_div(self.mesh, ampl_norm, sxy,
                          pngname=name + "_ampl_n.png")
        self.contourf_div(self.mesh, phas_norm, sxy,
                          pngname=name + "_phas_n.png")
        self.contourf_comp(self.mesh, ampl_norm, g_ampl, sxy,
                           pngname=name + "_ampl_compare.png")
        self.contourf_comp(self.mesh, db10, g_db10, sxy,
                           pngname=name + "_10db_compare.png", level=level)
        self.contourf_div(self.mesh, g_ampl, sxy,
                          pngname=name + "_ampl_g.png")
        self.contourf_div(self.mesh, g_db10, sxy,
                          pngname=name + "_10db_g.png", level=level)


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

    obj.mesh = rot_mesh(obj.mesh, rot=np.deg2rad(5.0))
    obj.setup_func()
    obj.contourf_div(obj.mesh, np.abs(obj.func), pngname=obj.tempname + "_rot")
