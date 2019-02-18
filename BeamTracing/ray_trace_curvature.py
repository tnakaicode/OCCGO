from OCCUtils.Topology import Topo
from OCCUtils.Construct import project_edge_onto_plane, project_point_on_curve
from OCCUtils.Construct import make_wire, make_edge, make_plane, make_line, make_loft
from OCCUtils.Construct import dir_to_vec, vec_to_dir
from OCC.BRep import BRep_Tool
from OCC.TopLoc import TopLoc_Location
from OCC.TopoDS import TopoDS_Face
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
from OCC.gp import gp_Mat
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
    p1 = gp_Pnt()
    vx, vy = gp_Vec(), gp_Vec()
    vxx, vyy = gp_Vec(), gp_Vec()
    vxy = gp_Vec()
    GeomLProp_SurfaceTool.D2(h_surf, u, v, p1, vx, vy, vxx, vyy, vxy)
    vx.Normalize()
    vy.Normalize()
    vxx.Normalize()
    vyy.Normalize()
    vxy.Normalize()
    print(p1)
    print("vx", vx)
    print("vy", vy)
    print("vxx", vxx)
    print("vyy", vyy)
    print("vxy", vxy)
    print(vx.Dot(vyy))
    print(vy.Dot(vxx))


if __name__ == "__main__":
    from src.RayTrace.RaySystem import RaySystem, SurfSystem, OptSystem, Multi_RaySystem
    from src.RayTrace.ray_setup import get_axs, get_deg
    from src.Unit import convert_SI, convert
    from src.geomtory import curvature
    from src.pyocc.surface import surf_spl

    argvs = sys.argv
    parser = OptionParser()
    parser.add_option("--dir", dest="dir", default="./")
    opt, argc = parser.parse_args(argvs)
    print(argc, opt)

    init = "surf1"
    surf = ["surf2", "surf3", "surf4"]

    surf1 = SurfSystem("./", "surf1")
    surf2 = SurfSystem("./", "surf2")

    h_surf = BRep_Tool.Surface(surf1.srf)
    second_derivative(h_surf, 0.5, 0.5)
    #second_derivative(h_surf, 1.0, 0)

    h_surf = BRep_Tool.Surface(surf2.srf)
    second_derivative(h_surf, 0.5, 0.5)
    #second_derivative(h_surf, 1.0, 0)

    print(surf1.name, surf1.axs.Location())
    print(surf2.name, surf2.axs.Location())
    print(surf2.srf)
    loc_surf = surf2.srf
    loc_trsf = gp_Trsf()
    loc_trsf.SetTransformation(gp_Ax3(), surf2.axs)
    loc_face = TopLoc_Location(loc_trsf)
    loc_surf.Move(loc_face)
    h_surf = BRep_Tool.Surface(loc_surf)
    second_derivative(h_surf, 0.5, 0.5)
    #second_derivative(h_surf, 1.0, 0)

    h_surf = BRep_Tool.Surface(surf2.srf)
    second_derivative(h_surf, 0.5, 0.5)
