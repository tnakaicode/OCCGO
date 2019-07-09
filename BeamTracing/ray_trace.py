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
from src.RayTrace.ray_setup import get_axs, get_deg
from src.Unit import convert_SI, convert

from OCCUtils.Topology import Topo
from OCCUtils.Construct import project_edge_onto_plane, project_point_on_curve
from OCCUtils.Construct import make_wire, make_edge, make_plane, make_line, make_loft
from OCCUtils.Construct import dir_to_vec, vec_to_dir
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

if __name__ == "__main__":
    argvs = sys.argv
    parser = OptionParser()
    parser.add_option("--dir", dest="dir", default="./")
    opt, argc = parser.parse_args(argvs)
    print(argc, opt)

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

    ax = obj.tar.Move_Axs(obj.tar.beam, obj.tar.axs, gp_Ax3())
    ax0 = obj.tar.Move_Axs(obj.tar.beam_rght, obj.tar.axs, gp_Ax3())
    ax1 = obj.tar.Move_Axs(obj.tar.beam_left, obj.tar.axs, gp_Ax3())
    ax2 = obj.tar.Move_Axs(obj.tar.beam_uppr, obj.tar.axs, gp_Ax3())
    ax3 = obj.tar.Move_Axs(obj.tar.beam_bott, obj.tar.axs, gp_Ax3())
    print(ax.Location())
    print(ax0.Location())
    print(ax1.Location())
    print(ax2.Location())
    print(ax3.Location())

    obj.Display_Shape(["BLUE", "GREEN"])

    print(obj.tar.beam.Location())

    for idx, name in enumerate(surf[:-1]):
        print(name)
        obj.ini = obj.tar
        obj.tar = SurfSystem("./", surf[idx + 1])
        print(obj.ini.beam.Location())

        obj.MultiReflect()
        print(obj.tar.beam.Location())
        print(obj.tar.beam_rght.Location())
        print(obj.tar.beam_left.Location())
        print(obj.tar.beam_uppr.Location())
        print(obj.tar.beam_bott.Location())

        obj.Display_Shape(["BLUE", "GREEN"])

    obj.display.FitAll()
    obj.start_display()
