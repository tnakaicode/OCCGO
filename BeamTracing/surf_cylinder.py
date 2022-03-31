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

from src.fileout import occ_to_grasp_cor, occ_to_grasp_rim
from src.geomtory import curvature, grasp_sfc
from src.geomtory import fit_surf
from src.pyocc.surface import surf_spl
from src.pyocc.export import export_STEPFile_single

from OCC.Display.SimpleGui import init_display
from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_Circ, gp_Circ2d
from OCC.Core.Geom import Geom_Plane, Geom_Surface, Geom_BSplineSurface
from OCC.Core.Geom import Geom_Curve, Geom_Line, Geom_Ellipse
from OCC.Core.Geom import Geom_Circle
from OCC.Core.BRep import BRep_Tool
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeFace
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeEdge
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeWire
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeShell
from OCC.Core.BRepOffsetAPI import BRepOffsetAPI_ThruSections
from OCC.Core.BRepOffsetAPI import BRepOffsetAPI_MakePipe
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakePrism
from OCCUtils.Construct import vec_to_dir, dir_to_vec
from OCCUtils.Construct import make_edge


def wxy_wire(pnt, wxy=[10, 20]):
    if wxy[0] >= wxy[1]:
        ax2 = gp_Ax2(pnt, gp_Dir(0, 0, 1), gp_Dir(1, 0, 0))
        w_x = wxy[0]
        w_y = wxy[1]
    elif wxy[1] >= wxy[0]:
        ax2 = gp_Ax2(pnt, gp_Dir(0, 0, 1), gp_Dir(0, 1, 0))
        w_x = wxy[1]
        w_y = wxy[0]
    else:
        ax2 = gp_Ax2(pnt, gp_Dir(0, 0, 1), gp_Dir(1, 0, 0))
        w_x = wxy[0]
        w_y = wxy[1]
    obj = Geom_Ellipse(ax2, w_x, w_y).Elips()
    obj = BRepBuilderAPI_MakeWire(BRepBuilderAPI_MakeEdge(obj).Edge()).Wire()
    return obj


if __name__ == "__main__":
    argvs = sys.argv
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", dest="dir", default="./")
    parser.add_argument("--surf", dest="surf", default="cylinder")
    parser.add_argument("--radi", dest="radi", default=(10, 10),
                      type="float", nargs=2)
    parser.add_argument("--lxy", dest="lxy", default=(0, 10),
                      type="float", nargs=2)
    parser.add_argument("--rxy", dest="rxy", default=(0, 0),
                      type="float", nargs=2)
    opt = parser.parse_args()
    print(argc, opt)

    display, start_display, add_menu, add_function_to_menu = init_display()

    api = BRepOffsetAPI_ThruSections()

    pt = np.linspace(*opt.lxy, 10)
    pr_x = np.tan(np.deg2rad(opt.rxy[0])) * pt + opt.radi[0]
    pr_y = np.tan(np.deg2rad(opt.rxy[1])) * pt + opt.radi[1]
    for i, d in enumerate(pt):
        pnt = gp_Pnt(0, 0, pt[i])
        d_z = gp_Dir(0, 0, 1)
        wxy = [pr_x[i], pr_y[i]]
        obj = wxy_wire(pnt, wxy)
        display.DisplayShape(obj)
        api.AddWire(obj)

    api.Build()
    surf = api.Shape()
    display.DisplayShape(surf)

    export_STEPFile_single(surf, opt.dir + opt.surf + ".stp")

    display.FitAll()
    start_display()
