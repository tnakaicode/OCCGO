from OCCUtils.Construct import vec_to_dir, dir_to_vec
from OCC.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.gp import gp_Pnt, gp_Vec, gp_Dir
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
from optparse import OptionParser
sys.path.append(os.path.join('../'))


if __name__ == "__main__":
    from src.RayTrace.ray_setup import get_axs, get_deg
    from src.fileout import occ_to_grasp_cor, occ_to_grasp_rim

    argvs = sys.argv
    parser = OptionParser()
    parser.add_option("--dir", dest="dir", default="./")
    parser.add_option("--out", dest="out", default="./")
    parser.add_option("--refe", dest="refe", default="surf1")
    parser.add_option("--surf", dest="surf", default="surf2")
    parser.add_option("--pxyz", dest="pxyz",
                      default=(0, 0, 0), type="float", nargs=3)
    parser.add_option("--rxyz", dest="rxyz",
                      default=(0, 0, 0), type="float", nargs=3)
    opt, argc = parser.parse_args(argvs)
    print(argc, opt)

    filename = opt.dir + opt.refe + ".cor"
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
    occ_to_grasp_cor(ax, opt.surf, opt.out + opt.surf + ".cor")
