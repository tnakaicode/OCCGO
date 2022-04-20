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

from src.RayTrace.RaySystem import RaySystem, SurfSystem, OptSystem, Multi_RaySystem
from src.RayTrace.ray_setup import get_axs, get_deg
from src.Unit import convert_SI, convert
from src.geomtory import curvature
from src.pyocc.surface import surf_spl

from OCCUtils.Topology import Topo
from OCCUtils.Construct import project_edge_onto_plane, project_point_on_curve
from OCCUtils.Construct import make_wire, make_edge, make_plane, make_line, make_loft
from OCCUtils.Construct import dir_to_vec, vec_to_dir
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


if __name__ == "__main__":
    argvs = sys.argv
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", dest="dir", default="./")
    opt = parser.parse_args()
    print(opt, argvs)

    init = "surf1"
    surf = ["surf2", "surf3", "surf4"]

    obj = Multi_RaySystem("./", init, "surf2")
    obj.ini.beam = get_axs("./" + obj.ini.name + "_beam.cor", obj.ini.axs)

    obj.ini.beam_rght = get_axs(
        "./" + obj.ini.name + "_beam_rght.cor", obj.ini.axs)
    obj.ini.beam_left = get_axs(
        "./" + obj.ini.name + "_beam_left.cor", obj.ini.axs)
    obj.ini.beam_uppr = get_axs(
        "./" + obj.ini.name + "_beam_uppr.cor", obj.ini.axs)
    obj.ini.beam_bott = get_axs(
        "./" + obj.ini.name + "_beam_bott.cor", obj.ini.axs)

    obj.MultiReflect()
    print(obj.tar.beam.Location())
    print(obj.tar.beam_rght.Location())
    print(obj.tar.beam_left.Location())
    print(obj.tar.beam_uppr.Location())
    print(obj.tar.beam_bott.Location())

    rxy = [-10, -100]

    px = np.linspace(-1, 1, 100) * 25
    py = np.linspace(-1, 1, 100) * 25
    mesh = np.meshgrid(px, py)

    rx_0 = curvature(mesh[0], rxy[0], 0)
    ry_0 = curvature(mesh[1], rxy[1], 0)
    ph_0 = wavefront_xyz(*mesh, rx_0 + ry_0, obj.ini.beam)

    dist = obj.ini.beam.Location().Distance(obj.tar.beam.Location())

    rx_1 = curvature(mesh[0], rxy[0] + dist / 2, 0)
    ry_1 = curvature(mesh[1], rxy[1] + dist / 2, 0)
    ph_1 = wavefront_xyz(*mesh, rx_1 + ry_1 + dist / 2, obj.ini.beam)

    rx_r = curvature(mesh[0], rxy[0] + dist, 0)
    ry_r = curvature(mesh[1], rxy[1] + dist, 0)
    ph_r = wavefront_xyz(*mesh, rx_r + ry_r + dist, obj.ini.beam)

    obj.Display_Shape(["BLUE", "GREEN"])
    obj.display.DisplayShape(ph_0, color="RED")
    obj.display.DisplayShape(ph_1, color="RED")
    obj.display.DisplayShape(ph_r, color="YELLOW")

    print(obj.tar.beam.Location())

    obj.display.FitAll()
    obj.start_display()
