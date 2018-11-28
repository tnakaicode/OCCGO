import numpy as np
import matplotlib.pyplot as plt
import scipy.special
import json, glob
import sys, time, os
import scipy.constants as cnt
from unwrap.unwrap import unwrap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from linecache import getline, clearcache
from optparse import OptionParser  
from scipy.integrate import simps
from scipy.optimize import curve_fit
from scipy.special import jv, jvp, jn_zeros, jnp_zeros
from scipy.special import yv, yvp, yn_zeros, ynp_zeros
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "\..")

from src.Unit import convert_SI, convert

if __name__ == "__main__":
    argvs = sys.argv  
    parser = OptionParser()
    parser.add_option("--freq", dest="freq", default="170GHz")
    parser.add_option("--radi", dest="radi", default=(10.0, 25.0), type="float", nargs=2)
    parser.add_option("--mn" , dest="mn" , default=(31,11), type="int", nargs=2)
    parser.add_option("--nxy", dest="nxy", default=(200,200), type="int", nargs=2)
    parser.add_option("--lxy", dest="lxy", default=(100,100), type="float")
    opt, argc = parser.parse_args(argvs)
    print (argc, opt)

    freq = convert_SI (opt.freq, unit_in='GHz', unit_out='Hz')
    wave = cnt.c / freq * convert(unit_in="m", unit_out="mm")

    px = np.linspace (-opt.lxy[0]/2, opt.lxy[0]/2, opt.nxy[0])
    py = np.linspace (-opt.lxy[1]/2, opt.lxy[1]/2, opt.nxy[1])
    mesh = np.meshgrid (px, py)
    x, y = mesh[0], mesh[1]
    r, t = np.sqrt(x**2+y**2), np.arctan(y/x)
    m, n = opt.mn
    a, b = opt.radi

    if opt.mn[0] == 0:
        e_m = 1
    else:
        e_m = 2

    x0_i = jn_zeros (*opt.mn)[-1]
    x1_i = jnp_zeros (*opt.mn)[-1]

    func = jv(m, x1_i*r/b)*yvp(m, x1_i) - jvp(m, x1_i)*yv(m, x1_i*r/b)

    plt.figure()
    plt.contourf (x, y, func)
    plt.show()
