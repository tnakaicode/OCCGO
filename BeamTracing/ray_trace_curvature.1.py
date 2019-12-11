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

from src.RayTrace.RaySystem import RaySystem, SurfSystem, OptSystem, Multi_RaySystem
from src.RayTrace.ray_setup import get_axs, get_deg, axs_pln
from src.Unit import convert_SI, convert
from src.geomtory import curvature
from src.pyocc.surface import surf_spl, curv_spl

from OCCUtils.Topology import Topo
from OCCUtils.Construct import project_edge_onto_plane, project_point_on_curve
from OCCUtils.Construct import make_wire, make_edge, make_plane, make_line, make_loft
from OCCUtils.Construct import dir_to_vec, vec_to_dir
from OCC.Core.BRep import BRep_Tool, BRep_PointsOnSurface
from OCC.Core.BRep import BRep_ListNodeOfListOfPointRepresentation
from OCC.Core.TopLoc import TopLoc_Location
from OCC.Core.TopoDS import TopoDS_Face
from OCC.Core.GeomLProp import GeomLProp_SurfaceTool, GeomLProp_SLProps
from OCC.Core.GeomAbs import GeomAbs_C2, GeomAbs_C0, GeomAbs_G1, GeomAbs_G2
from OCC.Core.GeomAPI import GeomAPI_ProjectPointOnSurf, GeomAPI_ProjectPointOnCurve
from OCC.Core.GeomAPI import GeomAPI_PointsToBSplineSurface, GeomAPI_IntCS
from OCC.Core.GeomAPI import GeomAPI_ExtremaCurveCurve, GeomAPI_Interpolate
from OCC.Core.TColgp import TColgp_Array1OfPnt, TColgp_Array2OfPnt
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeFace
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_ParametersOutOfRange
from OCC.Core.Geom import Geom_Curve, Geom_Line, Geom_Ellipse
from OCC.Core.Geom import Geom_Plane, Geom_Surface, Geom_BSplineSurface, Geom_BSplineCurve
from OCC.Core.Geom import Geom_ConicalSurface, Geom_Conic
from OCC.Core.gp import gp_Pln, gp_Trsf, gp_Lin, gp_Elips, gp_Elips2d
from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir, gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_Mat
from OCC.Display.SimpleGui import init_display


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
    pu, pv = gp_Vec(), gp_Vec()
    puu, pvv = gp_Vec(), gp_Vec()
    puv = gp_Vec()
    prop = GeomLProp_SLProps(h_surf, u, v, 1, 1)
    GeomLProp_SurfaceTool.D2(h_surf, u, v, p1, pu, pv, puu, pvv, puv)
    e0 = pu.Crossed(pv)
    pu.Normalize()
    pv.Normalize()
    e0.Normalize()
    puu.Normalize()
    pvv.Normalize()
    puv.Normalize()
    print(p1)
    print("pu", pu)
    print("pv", pv)
    print("e0", e0)
    print("puu", puu)
    print("pvv", pvv)
    print("puv", puv)

    first_form = np.array([
        [pu.Dot(pu), pu.Dot(pv)],
        [pv.Dot(pu), pv.Dot(pv)]
    ])
    secnd_form = np.array([
        [e0.Dot(puu), e0.Dot(puv)],
        [e0.Dot(puv), e0.Dot(pvv)]
    ])

    print(first_form)
    print(secnd_form)
    print(prop.GaussianCurvature())
    print(prop.MeanCurvature())
    d1, d2 = gp_Dir(), gp_Dir()
    prop.CurvatureDirections(d1, d2)
    a1 = gp_Ax3()
    v1 = dir_to_vec(d1)
    v2 = dir_to_vec(d2)
    if pu.IsParallel(v1, 1 / 1000):
        c1 = prop.MaxCurvature()
        c2 = prop.MinCurvature()
        print(v1.Dot(pu), v1.Dot(pv))
        print(v2.Dot(pu), v2.Dot(pv))
    else:
        c1 = prop.MinCurvature()
        c2 = prop.MaxCurvature()
        print(v1.Dot(pu), v1.Dot(pv))
        print(v2.Dot(pu), v2.Dot(pv))
    print(c1, 1 / c1)
    print(c2, 1 / c2)

    px = np.linspace(-1, 1, 100) * 100
    p1_y = px**2 / c1
    p2_y = px**2 / c1
    curv1 = curv_spl(px, p1_y)
    curv2 = curv_spl(px, p2_y)


if __name__ == "__main__":
    argvs = sys.argv
    parser = OptionParser()
    parser.add_option("--dir", dest="dir", default="./")
    opt, argc = parser.parse_args(argvs)
    print(argc, opt)

    display, start_display, add_menu, add_function_to_menu = init_display()

    axs = gp_Ax3()
    deg = 30
    rad = 0
    obj = Geom_ConicalSurface(axs, np.deg2rad(deg), rad)
    surf = BRepBuilderAPI_MakeFace(
        obj.GetHandle(), -np.pi / 2, np.pi / 2, -100, 100, 1e-6).Face()

    display.DisplayShape(surf)
    display.DisplayShape(obj)
    display.DisplayShape(axs_pln(axs))

    display.FitAll()
    start_display()

    init = "surf1"
    surf = ["surf2", "surf3", "surf4"]

    surf1 = SurfSystem("./", "surf1")
    surf2 = SurfSystem("./", "surf2")

    h_surf = BRep_Tool.Surface(surf1.srf)
    second_derivative(h_surf, 0.5, 0.5)
    second_derivative(h_surf, 0.5, 0.0)
    second_derivative(h_surf, 0.0, 0.5)

    """h_surf = BRep_Tool.Surface(surf2.srf)
    second_derivative(h_surf, 0.5, 0.5)
    second_derivative(h_surf, 0.5, 0.0)
    second_derivative(h_surf, 0.0, 0.5)
    second_derivative(h_surf, 0.0, 0.0)"""

    """print(surf1.name, surf1.axs.Location())
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
    second_derivative(h_surf, 0.5, 0.5)"""
