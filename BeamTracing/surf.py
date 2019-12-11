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

from src.fileout import occ_to_grasp_cor, occ_to_grasp_rim
from src.geomtory import curvature, grasp_sfc
from src.geomtory import fit_surf
from src.pyocc.surface import surf_spl
from src.pyocc.export import export_STEPFile_single

from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3

if __name__ == "__main__":
    argvs = sys.argv
    parser = OptionParser()
    parser.add_option("--dir", dest="dir", default="./")
    parser.add_option("--surf", dest="surf", default="surf1")
    parser.add_option("--lxy", dest="lxy",
                      default=(500, 500), type="float", nargs=2)
    parser.add_option("--nxy", dest="nxy",
                      default=(100.0, 100.0), type="int", nargs=2)
    parser.add_option("--sxy", dest="sxy", default=(0, 0),
                      type="float", nargs=2)
    parser.add_option("--rxy", dest="rxy",
                      default=(100, 200), type="float", nargs=2)
    opt, argc = parser.parse_args(argvs)
    print(argc, opt)

    px = np.linspace(-1, 1, opt.nxy[0]) * opt.lxy[0] / 2
    py = np.linspace(-1, 1, opt.nxy[1]) * opt.lxy[1] / 2
    mesh = np.meshgrid(px, py)
    curx = curvature(mesh[0], opt.rxy[0], opt.sxy[0])
    cury = curvature(mesh[1], opt.rxy[1], opt.sxy[1])
    surf = curx + cury

    grasp_sfc(mesh, surf, opt.dir + opt.surf + "_mat.sfc")
