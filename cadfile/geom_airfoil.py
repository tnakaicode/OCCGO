import numpy as np
import matplotlib.pyplot as plt
import json
import glob
import sys
import time
import os
import urllib
import scipy.constants as cnt
from optparse import OptionParser

from OCC.Display.SimpleGui import init_display
from OCC.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.gp import gp_Pln, gp_Trsf, gp_Lin


def get_airfol_data(filename="./dae51.dat"):
    fp = urllib.request.urlopen(filename)
    dat = fp.readlines()
    print(filename)
    print(dat[0].split()[0])
    print(dat[1].split())
    n0 = int(float(dat[1].split()[0]))
    n1 = int(float(dat[1].split()[1]))

    ns0 = 3
    ne0 = ns0+n0
    ns1 = 3+n0+1
    ne1 = ns1+n1
    upp = []
    bot = []
    for line in dat[ns0:ne0]:
        x, y = [float(v) for v in line.split()]
        upp.append(np.array([x, y, 0]))

    for line in dat[ns1:ne1]:
        x, y = [float(v) for v in line.split()]
        bot.append(np.array([x, y, 0]))
    return np.array(upp), np.array(bot)


if __name__ == "__main__":
    argvs = sys.argv
    parser = OptionParser()
    parser.add_option("--name", dest="name", default="dae51")
    opt, argc = parser.parse_args(argvs)
    print(argc, opt)

    cfg = json.load(open("airfoil.json", "r"))
    airfilo_url = cfg["url"]
    airfilo_dat = cfg[opt.name]
    filename = airfilo_url + airfilo_dat
    print(filename)

    # get_airfol_data()
    upp, bot = get_airfol_data(filename)
    # upp.shape()
    print(upp[:, 0])

    plt.figure()
    plt.plot(upp[:, 0], upp[:, 1])
    plt.plot(bot[:, 0], bot[:, 1])
    plt.show()
