import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
from linecache import getline, clearcache
from optparse import OptionParser

sys.path.append(os.path.join("../"))
from src.base import plotocc, spl_face

import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)

from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.TopoDS import TopoDS_Shape, TopoDS_Compound
from OCC.Core.TopLoc import TopLoc_Location
from OCC.Core.BRep import BRep_Tool
from OCC.Core.Geom import Geom_Line
from OCC.Core.GeomLProp import GeomLProp_SLProps
from OCC.Core.GeomAPI import GeomAPI_ProjectPointOnSurf, GeomAPI_ProjectPointOnCurve
from OCC.Core.GeomAPI import GeomAPI_IntCS, GeomAPI_IntSS
from OCCUtils.Topology import Topo
from OCCUtils.Construct import make_box, make_line, make_wire, make_edge
from OCCUtils.Construct import make_plane, make_polygon
from OCCUtils.Construct import point_to_vector, vector_to_point
from OCCUtils.Construct import dir_to_vec, vec_to_dir

if __name__ == '__main__':
    argvs = sys.argv
    parser = OptionParser()
    parser.add_option("--dir", dest="dir", default="./")
    parser.add_option("--pxyz", dest="pxyz",
                      default=[0.0, 0.0, 0.0], type="float", nargs=3)
    opt, argc = parser.parse_args(argvs)
    print(opt, argc)

    obj = plotocc(touch=True)

    ax0 = gp_Ax3(gp_Pnt(100, 100, 20),
                 gp_Dir(0.1, 0.2, 1.0),
                 gp_Dir(0.0, 0.5, 1.0))
    px = np.linspace(-1, 1, 500) * 100 + 50
    py = np.linspace(-1, 1, 700) * 100 - 50
    mesh = np.meshgrid(px, py)
    data = mesh[0]**2 / 1000 + mesh[1]**2 / 2000
    surf = spl_face(*mesh, data, ax0)

    ax1 = gp_Ax3(gp_Pnt(150, 150, -20),
                 gp_Dir(0.5, 0.1, 1.0),
                 gp_Dir(0.0, 0.1, 1.0))
    pnt = ax1.Location()
    api = GeomAPI_ProjectPointOnSurf(pnt, BRep_Tool.Surface(surf))
    print(api.NbPoints(), api.NearestPoint())
    for i in range(api.NbPoints()):
        idx = i + 1
        u, v = api.Parameters(idx)
        obj.display.DisplayShape(api.Point(idx))
        print(idx, u, v)

    print(GeomAPI_IntCS(Geom_Line(ax1.Axis()), BRep_Tool.Surface(surf)).NbPoints())
    u, v, w = GeomAPI_IntCS(Geom_Line(ax1.Axis()), BRep_Tool.Surface(surf)).Parameters(1)
    api = BRep_Tool.Surface(surf)
    p, vx, vy = gp_Pnt(), gp_Vec(), gp_Vec()
    #print(api.D1(u, v))
    
    api = GeomLProp_SLProps(BRep_Tool.Surface(surf), u, v, 1, 0.1E-03)
    print(api.GaussianCurvature())
    print(api.MinCurvature(), api.MeanCurvature(), api.MaxCurvature())
    print(dir_to_vec(api.Normal()))
    du, dv = gp_Dir(), gp_Dir()
    api.TangentU(du)
    api.TangentV(dv)
    print(dir_to_vec(du))
    print(dir_to_vec(dv))
    dx, dy = gp_Dir(), gp_Dir()
    api.CurvatureDirections(dx, dy)
    print(dir_to_vec(dx))
    print(dir_to_vec(dy))
    
    obj.show_axs_pln(ax0, scale=100, name="axis-0")
    obj.show_axs_pln(ax1, scale=100, name="axis-1")
    obj.display.DisplayShape(surf, color="BLUE", transparency=0.9)
    obj.show_axs_pln(scale=100)
    obj.show()
