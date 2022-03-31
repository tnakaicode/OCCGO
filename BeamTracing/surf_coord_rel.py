import numpy as np
import matplotlib.pyplot as plt
import json
import glob
import sys
import time
import os
from unwrap.unwrap import unwrap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from linecache import getline, clearcache
from scipy.integrate import simps
import argparse
sys.path.append(os.path.join('../'))

from src.RayTrace.ray_setup import get_axs, get_deg
from src.fileout import occ_to_grasp_cor, occ_to_grasp_rim

from OCCUtils.Construct import vec_to_dir, dir_to_vec
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir

if __name__ == "__main__":
    argvs = sys.argv
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", dest="dir", default="./")
    parser.add_argument("--surf", dest="surf", default="surf1")
    parser.add_argument("--pxyz", dest="pxyz",
                      default=(0, 0, 0), type=float, nargs=3)
    parser.add_argument("--rxyz", dest="rxyz",
                      default=(0, 0, 0), type=float, nargs=3)
    opt = parser.parse_args()
    print(argc, opt)

    filename = opt.dir + opt.surf + ".cor"
    ax = get_axs(filename)
    print(dir_to_vec(ax.Direction()), ax.Location())
    ax.Translate(gp_Vec(gp_Pnt(), gp_Pnt(*opt.pxyz)))
    for i, deg in enumerate(opt.rxyz):
        if i == 0:
            axs = gp_Ax1(ax.Location(), ax.XDirection())
        elif i == 1:
            axs = gp_Ax1(ax.Location(), ax.YDirection())
        elif i == 2:
            axs = gp_Ax1(ax.Location(), ax.Direction())
        else:
            axs = gp_Ax1(ax.Location(), ax.Direction())
        ax.Rotate(axs, np.deg2rad(deg))
    print(dir_to_vec(ax.Direction()), ax.Location())
    occ_to_grasp_cor(ax, opt.surf, filename)
