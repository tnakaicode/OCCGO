import shutil
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
from scipy.stats import multivariate_normal
from linecache import getline, clearcache
from optparse import OptionParser

from sympy import im

sys.path.append(os.path.join("../"))
from src.base import plot2d, plot3d
from src.Unit import knum_from_freq, knum_from_wave
from src.profile import moment, gcf_calc, normalize_integrate
from src.profile import gaussian_func, rot_mesh, get_cov, get_centroid
from src.geomtory import curvature
from src.Gaussian import GaussianProfile, gen_noise, copy_file

import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)


def moment_ensable(mesh, func, g_func):
    func1 = normalize_integrate(mesh, func)
    func2 = normalize_integrate(mesh, g_func)
    #func2 = g_func
    fg_func = normalize_integrate(mesh, func1 * func2)
    #fg_func = func * g_func
    sx = moment(mesh, fg_func, [1, 0])
    sy = moment(mesh, fg_func, [0, 1])
    wx = moment(mesh, fg_func, [2, 0]) * 2
    rt = moment(mesh, fg_func, [1, 1]) * 2
    wy = moment(mesh, fg_func, [0, 2]) * 2
    sxy = np.array([sx, sy])
    cov = np.array([
        [wx, rt],
        [rt, wy]
    ])
    return sxy, cov


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

    ampl = obj.ampl
    mesh_xy = np.stack(obj.mesh, -1)
    sxy, cov = moment_ensable(obj.mesh, ampl, ampl)
    g_ampl = multivariate_normal.pdf(mesh_xy, mean=sxy, cov=cov)
    wx, wy = np.sqrt(cov[0, 0]), np.sqrt(cov[1, 1])
    rot = np.rad2deg(np.arccos(cov[0, 1] / wx / wy))
    gcf = gcf_calc(obj.mesh, ampl, g_ampl)
    print(sxy, wx, wy, rot, gcf)
    for i in range(100):
        g_ampl = multivariate_normal.pdf(mesh_xy, mean=sxy, cov=cov)
        sxy, cov = moment_ensable(obj.mesh, ampl, g_ampl)
        wx, wy = np.sqrt(cov[0, 0]), np.sqrt(cov[1, 1])
        rot = np.rad2deg(np.arccos(cov[0, 1] / wx / wy))
        gcf = gcf_calc(obj.mesh, ampl, g_ampl)
        print(sxy, wx, wy, rot, gcf)
        wxy, mat = np.linalg.eig(cov)
        rot = np.rad2deg(np.arctan2(mat[0, 1], mat[1, 1]))
        print(np.sqrt(wxy), rot)
