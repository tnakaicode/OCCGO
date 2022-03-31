import numpy as np
import matplotlib.pyplot as plt
import json
import sys
import time
import os
import glob
from unwrap.unwrap import unwrap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from linecache import getline, clearcache
from scipy.integrate import simps
import argparse
sys.path.append(os.path.join('../'))

from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_Pln, gp_Trsf, gp_Lin
from OCC.Core.Geom import Geom_Plane, Geom_Surface, Geom_BSplineSurface
from OCC.Core.Geom import Geom_Curve, Geom_Line
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeFace
from OCC.Core.TColgp import TColgp_Array1OfPnt, TColgp_Array2OfPnt
from OCC.Core.GeomAPI import GeomAPI_PointsToBSplineSurface, GeomAPI_IntCS
from OCC.Core.GeomAPI import GeomAPI_ProjectPointOnSurf, GeomAPI_ProjectPointOnCurve
from OCC.Core.GeomAbs import GeomAbs_C2, GeomAbs_C0, GeomAbs_G1, GeomAbs_G2
from OCC.Core.GeomLProp import GeomLProp_SurfaceTool
from OCC.Core.TopLoc import TopLoc_Location
from OCC.Core.BRep import BRep_Tool
from OCCUtils.Construct import dir_to_vec, vec_to_dir
from OCCUtils.Construct import make_wire, make_edge, make_plane, make_line
from OCCUtils.Construct import project_edge_onto_plane, project_point_on_curve
from OCCUtils.Topology import Topo
from OCC.Display.SimpleGui import init_display

from src.pyocc.load import read_step_file_shapes, read_step_file, read_step, read_iges
from src.pyocc.export import ExportCAFMethod, ExportMethod
from src.pyocc.OCCQt import Viewer
from src.RayTrace.ray_setup import get_axs, get_deg, axs_pln, reflect


def Transform(ax0=gp_Ax3(), ax1=gp_Ax3(), shape=[]):
    trf = gp_Trsf()
    trf.SetTransformation(ax1, ax0)
    loc_face = TopLoc_Location(trf)

    for shp in shape:
        shp.Location(loc_face)


def SetViewer():
    viewer = Viewer()
    viewer.on_select()


if __name__ == "__main__":
    argvs = sys.argv
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", dest="dir", default="../BeamTracing/")
    parser.add_argument("--out", dest="out", default="./")
    parser.add_argument("--file", dest="file", default="surf1.stp")
    parser.add_argument("--pxyz", dest="pxyz",
                      default=(0, 0, 0), type="float", nargs=3)
    parser.add_argument("--rxyz", dest="rxyz",
                      default=(0, 0, 0), type="float", nargs=3)
    parser.add_argument("--fileout", dest="fileout", default="surf2.stp")
    opt = parser.parse_args()
    print(argc, opt)

    idx = opt.file.split(".")[-1]

    if idx in ["stp", "step"]:
        shpe = read_step_file_shapes(opt.dir + opt.file)
    elif idx in ["igs", "iges"]:
        shpe = read_iges(opt.dir + opt.file)
    else:
        print("Incorrect file index")
        sys.exit(0)

    print(shpe)
    ax = gp_Ax3(gp_Pnt(*opt.pxyz), gp_Dir(0, 0, 1))
    for i, deg in enumerate(opt.rxyz):
        if i == 0:
            axs = gp_Ax1(ax.Location(), ax.XDirection())
        elif i == 1:
            axs = gp_Ax1(ax.Location(), ax.YDirection())
        elif i == 2:
            axs = gp_Ax1(ax.Location(), ax.Direction())
        else:
            axs = gp_Ax1(ax.Location(), ax.Direction())
        ax.Rotate(axs, np.deg2rad(deg))
    Transform(gp_Ax3(), ax, shpe)

    display, start_display, add_menu, add_function_to_menu = init_display()

    viewer = Viewer()
    viewer.on_select()

    #add_menu ("Viewer")
    #add_function_to_menu ("Viewer", SetViewer)

    ex_obj = ExportMethod()
    for shp in shpe:
        ex_obj.add_shpe(shp)
    ex_obj.fileout(opt.out + opt.fileout)

    SetViewer()
    display.DisplayShape(shpe)
    display.DisplayShape(axs_pln(gp_Ax3()))
    display.DisplayShape(axs_pln(ax))
    display.DisplayShape(gp_Pnt())

    display.FitAll()
    start_display()
