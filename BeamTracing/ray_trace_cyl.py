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

from src.RayTrace.RaySystem import RaySystem, SurfSystem, OptSystem, Multi_RaySystem, set_surface
from src.RayTrace.ray_setup import get_axs, get_deg, axs_pln, reflect
from src.Unit import convert_SI, convert
from src.pyocc.export import export_STEPFile_single
from src.pyocc.load import read_step_file

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


def reflect_axs2(beam, surf, axs=gp_Ax3(), indx=1):
    p0, v0 = beam.Location(), dir_to_vec(beam.Direction())
    h_surf = BRep_Tool.Surface(surf)
    ray = Geom_Line(gp_Lin(p0, vec_to_dir(v0)))
    if GeomAPI_IntCS(ray.GetHandle(), h_surf).NbPoints() == 0:
        return beam, beam, None
    elif GeomAPI_IntCS(ray.GetHandle(), h_surf).NbPoints() == 1:
        return beam, beam, None
    GeomAPI_IntCS(ray.GetHandle(), h_surf).IsDone()
    u, v, w = GeomAPI_IntCS(ray.GetHandle(), h_surf).Parameters(indx)
    p1, vx, vy = gp_Pnt(), gp_Vec(), gp_Vec()
    GeomLProp_SurfaceTool.D1(h_surf, u, v, p1, vx, vy)
    vz = vx.Crossed(vy)
    vx.Normalize()
    vy.Normalize()
    vz.Normalize()
    v1 = v0.Mirrored(gp_Ax2(p1, vec_to_dir(vz), vec_to_dir(vx)))
    norm_ax = gp_Ax3(p1, vec_to_dir(vz), vec_to_dir(vx))
    beam_ax = gp_Ax3(p1, vec_to_dir(v1), beam.XDirection().Reversed())
    return beam_ax, norm_ax, 1


if __name__ == "__main__":
    argvs = sys.argv
    parser = OptionParser()
    parser.add_option("--dir", dest="dir", default="./")
    opt, argc = parser.parse_args(argvs)
    print(argc, opt)

    display, start_display, add_menu, add_function_to_menu = init_display()

    beam = get_axs("./beam_cyl.cor")
    surf = set_surface("./cylinder.stp")
    display.DisplayShape(surf)
    display.DisplayShape(axs_pln(gp_Ax3()))
    display.DisplayShape(axs_pln(beam))

    axs1 = beam
    val = 1
    while val != None:
        axs0 = axs1
        axs1, axs, val = reflect_axs2(axs0, surf, indx=2)
        get_deg(axs, dir_to_vec(axs1.Direction()))
        if val != None:
            display.DisplayShape(
                make_line(axs0.Location(), axs1.Location()), color="GREEN")

    display.DisplayShape(axs_pln(axs1))

    display.FitAll()
    start_display()
