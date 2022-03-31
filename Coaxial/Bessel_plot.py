import numpy as np
import matplotlib.pyplot as plt
import scipy.special
import json
import glob
import sys
import time
import os
import scipy.constants as cnt
from unwrap.unwrap import unwrap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from linecache import getline, clearcache
import argparse
from scipy.integrate import simps
from scipy.optimize import curve_fit
from scipy.special import jv, jvp, jn_zeros, jnp_zeros
from scipy.special import yv, yvp, yn_zeros, ynp_zeros

sys.path.append(os.path.join('..'))
from src.base import plot2d
from src.Unit import convert_SI, convert


def circle_wg_tm_mn(mesh, m=31, n=11, phi=90.0, radi=25.0):
    x, y = mesh[0], mesh[1]
    r, t = np.sqrt(x**2 + y**2), np.arctan(y / x)
    b = radi

    x0_mn = jn_zeros(m, n)[-1]
    x1_mn = jnp_zeros(m, n)[-1]

    if m == 0:
        e_m = 1
    else:
        e_m = 2

    func = np.sqrt(e_m / np.pi)
    func *= 1 / np.abs(jvp(m + 1, x0_mn))
    func *= (-jvp(m, x0_mn * r / b, 0) * np.cos(m * t) + m /
             x0_mn * jvp(m, x0_mn * r / b, 0) / r * np.sin(m * t))
    return func


def circle_wg_te_mn(mesh, m=31, n=11, phi=90.0, radi=25.0):
    x, y = mesh[0], mesh[1]
    r, t = np.sqrt(x**2 + y**2), np.arctan(y / x)
    b = radi

    x0_mn = jn_zeros(m, n)[-1]
    x1_mn = jnp_zeros(m, n)[-1]

    if m == 0:
        e_m = 1
    else:
        e_m = 2

    func = np.sqrt(e_m / np.pi)
    func *= 1 / (np.sqrt(x1_mn**2 - m**2) * np.abs(jvp(m + 1, x0_mn)))
    func *= (-m * jvp(m, x1_mn * r / b, 0) / r * np.cos(m * t) +
             x1_mn * jvp(m, x1_mn * r / b, 1) / b * np.sin(m * t))
    return func


def coaxial_wg_tm_mn(mesh, m=31, n=11, phi=90.0, radi_a=0.0, radi_b=25.0):
    x, y = mesh[0], mesh[1]
    r, t = np.sqrt(x**2 + y**2), np.arctan(y / x)
    a, b = radi_a, radi_b

    x0_mn = jn_zeros(m, n)[-1]
    x1_mn = jnp_zeros(m, n)[-1]

    if m == 0:
        e_m = 1
    else:
        e_m = 2

    func = np.sqrt(e_m / np.pi)
    func *= 1 / np.abs(jvp(m + 1, x0_mn))
    func *= (-jvp(m, x0_mn * r / b, 0) * np.cos(m * t) + m /
             x0_mn * jvp(m, x0_mn * r / b, 0) / r * np.sin(m * t))
    return func


if __name__ == "__main__":
    argvs = sys.argv
    parser = argparse.ArgumentParser()
    parser.add_argument("--freq", dest="freq", default="170GHz")
    parser.add_argument("--radi", dest="radi",
                      default=(10.0, 25.0), type=float, nargs=2)
    parser.add_argument("--mn", dest="mn", default=(31, 11), type=int, nargs=2)
    parser.add_argument("--nxy", dest="nxy",
                      default=(200, 200), type=int, nargs=2)
    parser.add_argument("--lxy", dest="lxy", default=25, type=float)
    opt = parser.parse_args()
    print(argc, opt)

    freq = convert_SI(opt.freq, unit_in='GHz', unit_out='Hz')
    wave = cnt.c / freq * convert(unit_in="m", unit_out="mm")

    px = np.linspace(-1, 1, opt.nxy[0]) * opt.lxy / 2
    py = np.linspace(-1, 1, opt.nxy[1]) * opt.lxy / 2
    mesh = np.meshgrid(px, py)
    x, y = mesh[0], mesh[1]
    r, t = np.sqrt(x**2 + y**2), np.arctan(y / x)
    m, n = opt.mn
    a, b = opt.radi

    if m == 0:
        e_m = 1
    else:
        e_m = 2

    x0_i = jn_zeros(m, n)[-1]
    x1_i = jnp_zeros(m, n)[-1]
    print(x0_i, x1_i)

    obj = plot2d(aspect="auto")
    obj.axs.contourf(x, y, circle_wg_tm_mn(mesh))
    obj.SavePng_Serial()

    obj.new_2Dfig(aspect="auto")
    obj.axs.contourf(x, y, circle_wg_te_mn(mesh))
    obj.SavePng_Serial()
