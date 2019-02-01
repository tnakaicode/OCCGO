from OCCUtils.Topology import Topo
from OCCUtils.Construct import project_edge_onto_plane, project_point_on_curve
from OCCUtils.Construct import make_wire, make_edge, make_plane, make_line, make_loft
from OCCUtils.Construct import dir_to_vec, vec_to_dir
from OCCUtils.Topology import Topo
from OCC.BRep import BRep_Tool
from OCC.TopLoc import TopLoc_Location
from OCC.GeomLProp import GeomLProp_SurfaceTool
from OCC.GeomAbs import GeomAbs_C2, GeomAbs_C0, GeomAbs_G1, GeomAbs_G2
from OCC.GeomAPI import GeomAPI_ProjectPointOnSurf, GeomAPI_ProjectPointOnCurve
from OCC.GeomAPI import GeomAPI_PointsToBSplineSurface, GeomAPI_IntCS
from OCC.TColgp import TColgp_Array1OfPnt, TColgp_Array2OfPnt
from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeFace
from OCC.Geom import Geom_Curve, Geom_Line, Geom_Ellipse
from OCC.Geom import Geom_Plane, Geom_Surface, Geom_BSplineSurface
from OCC.gp import gp_Pln, gp_Trsf, gp_Lin, gp_Elips, gp_Elips2d
from OCC.gp import gp_Pnt, gp_Vec, gp_Dir, gp_Ax1, gp_Ax2, gp_Ax3
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
from optparse import OptionParser
sys.path.append(os.path.join('../'))


if __name__ == "__main__":
    from src.plot import plot_contour_sub
    from src.geomtory import curvature

    argvs = sys.argv
    parser = OptionParser()
    parser.add_option("--dir", dest="dir", default="./")
    parser.add_option("--radi", dest="radi", default=(-1000, -1000),
                      type="float", nargs=2)
    opt, argc = parser.parse_args(argvs)
    print(argc, opt)

    dist = 100.0

    px = np.linspace(-1, 1, 100) * 100
    py = np.linspace(-1, 1, 100) * 100
    mesh = np.meshgrid(px, py)
    ph_0 = curvature(mesh[0], opt.radi[0], 0) + curvature(mesh[1], opt.radi[1], 0)
    ph_r = ph_0 + dist

    plot_contour_sub(mesh, ph_0, dirname="phas_r0")
    plot_contour_sub(mesh, ph_r, dirname="phas_r")
