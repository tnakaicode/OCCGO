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

from src.RayTrace.RaySystem import GOSystem
from src.RayTrace.ray_setup import get_axs, get_deg
from src.RayTrace.SurfSystem import GaussSystem
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


if __name__ == "__main__":
    argvs = sys.argv
    parser = OptionParser()
    parser.add_option("--dir", dest="dir", default="./")
    parser.add_option("--freq", dest="freq", default="170GHz")
    opt, argc = parser.parse_args(argvs)
    print(argc, opt)

    freq = convert_SI(opt.freq, unit_in='GHz', unit_out='Hz')
    wave = cnt.c / freq * convert(unit_in="m", unit_out="mm")

    obj = GOSystem("./", "surf1", "surf2", wave)
    obj.ini.Init_Beam()
    obj.tar.beam = obj.tar.axs

    obj.Reflect()
    obj.Display_Shape()
    obj.display.DisplayShape(obj.ini.wave)

    obj.ini = obj.tar
    obj.tar = GaussSystem ("./", "surf3")
    obj.Reflect()
    obj.Display_Shape()
    
    obj.display.FitAll()
    obj.start_display()
