from OCCUtils.Topology import Topo
from OCCUtils.Construct import project_edge_onto_plane, project_point_on_curve
from OCCUtils.Construct import make_wire, make_edge, make_plane, make_line, make_loft
from OCCUtils.Construct import dir_to_vec, vec_to_dir
from OCCUtils.Topology import Topo
from OCC.Core.BRep import BRep_Tool
from OCC.Core.TopLoc import TopLoc_Location
from OCC.Core.GeomLProp import GeomLProp_SurfaceTool
from OCC.Core.GeomAbs import GeomAbs_C2, GeomAbs_C0, GeomAbs_G1, GeomAbs_G2
from OCC.Core.GeomAPI import GeomAPI_ProjectPointOnSurf, GeomAPI_ProjectPointOnCurve
from OCC.Core.GeomAPI import GeomAPI_PointsToBSplineSurface, GeomAPI_IntCS
from OCC.Core.TColgp import TColgp_Array1OfPnt, TColgp_Array2OfPnt
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeFace
from OCC.Core.Geom import Geom_Curve, Geom_Line, Geom_Ellipse
from OCC.Core.Geom import Geom_Plane, Geom_Surface, Geom_BSplineSurface
from OCC.Core.gp import gp_Pln, gp_Trsf, gp_Lin, gp_Elips, gp_Elips2d
from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir, gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Display.SimpleGui import init_display

import numpy as np
import matplotlib.pyplot as plt
import json
import glob
import sys
import time
import os
import scipy.constants as cnt
from unwrap.unwrap import unwrap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from linecache import getline, clearcache
from scipy.integrate import simps
import argparse
sys.path.append(os.path.join('../'))


if __name__ == "__main__":
    from src.plot import plot_contour_sub
    from src.geomtory import curvature

    argvs = sys.argv
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", dest="dir", default="./")
    parser.add_argument("--radi", dest="radi", default=(-1000, -500),
                      type="float", nargs=2)
    parser.add_argument("--dist", dest="dist", default=500.0, type="float")
    opt = parser.parse_args()
    print(argc, opt)

    px = np.linspace(-1, 1, 100) * 100
    py = np.linspace(-1, 1, 100) * 100
    mesh = np.meshgrid(px, py)
    
    rx_0 = curvature(mesh[0], opt.radi[0], 0)
    ry_0 = curvature(mesh[1], opt.radi[1], 0)   
    ph_0 = rx_0 + ry_0
    
    rx_r = curvature(mesh[0], opt.radi[0] + opt.dist, 0)
    ry_r = curvature(mesh[1], opt.radi[1] + opt.dist, 0) 
    ph_r = ph_0 + rx_r + ry_r
    
    plot_contour_sub(mesh, ph_0, dirname="phas_r0")
    plot_contour_sub(mesh, ph_r, dirname="phas_r")
