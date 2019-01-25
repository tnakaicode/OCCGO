from OCC.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.gp import gp_Ax1, gp_Ax2, gp_Ax3
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
    from src.fileout import occ_to_grasp_cor, occ_to_grasp_rim
    from src.geomtory import curvature, grasp_sfc
    from src.geomtory import fit_surf
    from src.pyocc.surface import surf_spl
    from src.pyocc.export import export_STEPFile_single

    argvs = sys.argv
    parser = OptionParser()
    parser.add_option("--dir", dest="dir", default="./")
    parser.add_option("--surf", dest="surf", default="surf1")
    opt, argc = parser.parse_args(argvs)
    print(argc, opt)

    filename = opt.dir + opt.surf + "_mat.sfc"
    xs, ys, xe, ye = [float(v) for v in getline(filename, 2).split()]
    nx, ny = [int(v) for v in getline(filename, 3).split()]
    px = np.linspace (xs, xe, nx)
    py = np.linspace (ys, ye, ny)
    mesh = np.meshgrid(px, py)
    data = np.loadtxt(filename, skiprows=3)
    surf = surf_spl (*mesh, data)

    export_STEPFile_single(surf, opt.dir + opt.surf + ".stp")
