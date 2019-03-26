import numpy as np
import matplotlib.pyplot as plt
import json, glob
import sys, time, os
import scipy.constants as cnt
from unwrap.unwrap import unwrap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from linecache import getline, clearcache
from scipy.integrate import simps
from optparse import OptionParser
sys.path.append(os.path.join('../'))

from src.RayTrace.RaySystem import GaussSystem, GOSystem
from src.RayTrace.ray_setup import get_axs, get_deg
from src.Unit import convert_SI, convert
from src.geomtory import curvature
from src.pyocc.surface import surf_spl

from OCC.Display.SimpleGui import init_display
from OCC.gp import gp_Pln, gp_Trsf, gp_Lin, gp_Elips, gp_Elips2d
from OCC.gp import gp_Pnt, gp_Vec, gp_Dir, gp_Ax1, gp_Ax2, gp_Ax3
from OCC.gp import gp_Mat
from OCC.BRep import BRep_Tool, BRep_PointsOnSurface
from OCC.TopLoc import TopLoc_Location
from OCC.TopoDS import TopoDS_Face
from OCC.GeomLProp import GeomLProp_SurfaceTool, GeomLProp_SLProps
from OCC.GeomAbs import GeomAbs_C2, GeomAbs_C0, GeomAbs_G1, GeomAbs_G2
from OCC.GeomAPI import GeomAPI_ProjectPointOnSurf, GeomAPI_ProjectPointOnCurve
from OCC.GeomAPI import GeomAPI_PointsToBSplineSurface, GeomAPI_IntCS
from OCC.TColgp import TColgp_Array1OfPnt, TColgp_Array2OfPnt
from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeFace
from OCC.Geom import Geom_Curve, Geom_Line, Geom_Ellipse
from OCC.Geom import Geom_Plane, Geom_Surface, Geom_BSplineSurface
from OCC.Geom import Geom_ConicalSurface, Geom_Conic
from OCCUtils.Topology import Topo
from OCCUtils.Construct import project_edge_onto_plane, project_point_on_curve
from OCCUtils.Construct import make_wire, make_edge, make_plane, make_line, make_loft
from OCCUtils.Construct import dir_to_vec, vec_to_dir


def wavefront(rxy=[1000, 1000], axs=gp_Ax3()):
    px = np.linspace(-1, 1, 100) * 10
    py = np.linspace(-1, 1, 100) * 10
    mesh = np.meshgrid(px, py)

    rx_0 = curvature(mesh[0], rxy[0], 0)
    ry_0 = curvature(mesh[1], rxy[1], 0)
    ph_0 = rx_0 + ry_0
    phas = surf_spl(*mesh, ph_0)

    trf = gp_Trsf()
    trf.SetTransformation(axs, gp_Ax3())
    loc_face = TopLoc_Location(trf)
    phas.Location(loc_face)
    return phas


def wavefront_xyz(x, y, z, axs=gp_Ax3()):
    phas = surf_spl(x, y, z)

    trf = gp_Trsf()
    trf.SetTransformation(axs, gp_Ax3())
    loc_face = TopLoc_Location(trf)
    phas.Location(loc_face)
    return phas


def second_derivative(h_surf, u=0, v=0):
    prop = GeomLProp_SLProps(h_surf, u, v, 1, 1)

    d1, d2 = gp_Dir(), gp_Dir()
    prop.CurvatureDirections(d1, d2)
    v1 = dir_to_vec(d1)
    v2 = dir_to_vec(d2)
    c1 = prop.MaxCurvature()
    c2 = prop.MinCurvature()
    print("Max", c1, 1 / c1, v1)
    print("Min", c2, 1 / c2, v2)
    print(v1.Dot(v2))
    print(prop.Value())


if __name__ == "__main__":
    argvs = sys.argv
    parser = OptionParser()
    parser.add_option("--dir", dest="dir", default="./")
    opt, argc = parser.parse_args(argvs)
    print(argc, opt)

    obj = GOSystem("./", "surf1", "surf2")
    obj.ini.Init_Beam()

    h_surf = BRep_Tool.Surface(obj.ini.wave)
    second_derivative(h_surf, 0.5, 0.5)
