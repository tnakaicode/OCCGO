import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
from scipy.fftpack import fft, fft2
from scipy.fftpack import ifft, ifft2
from scipy.fftpack import fftfreq, fftshift
from unwrap.unwrap import unwrap
from optparse import OptionParser

sys.path.append(os.path.join("../"))
from base import plot2d

import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)

if __name__ == '__main__':
    argvs = sys.argv
    parser = OptionParser()
    parser.add_option("--dir", dest="dir", default="./")
    parser.add_option("--pxyz", dest="pxyz",
                      default=[0.0, 0.0, 0.0], type="float", nargs=3)
    opt, argc = parser.parse_args(argvs)
    print(opt, argc)

    obj = plot2d()

    px = np.linspace(-1, 1, 100) * 0.5 * np.pi + 0.1
    py = np.linspace(-1, 1, 200) * 0.5 * np.pi - 0.1
    mesh = np.meshgrid(px, py)
    #surf = mesh[0]**2 / 400 + mesh[0] / 20 + mesh[1] / 1000
    #surf *= mesh[0]**2 / 2000
    #surf = np.arctan2(mesh[0], mesh[1]) - 2 * np.arctan2(mesh[0], 2 * mesh[1])
    func = np.sin(np.exp(mesh[0]**2 / 10) * np.exp(mesh[1]**2 / 10 - 5))

    fft_func = fft2(func)
    fft_func *= np.exp(-10.0 * mesh[0]**2)
    fft_func *= np.exp(-15.0 * mesh[1]**2)
    g_func = ifft2(fft_func)
    phas = np.angle(g_func)

    ax1 = obj.add_axs(1, 2, 1, aspect="equal")
    ax2 = obj.add_axs(1, 2, 2, aspect="equal")
    ax1.contourf(*mesh, phas, cmap="jet")
    ax2.contourf(*mesh, unwrap(phas), cmap="jet")
    obj.SavePng()
