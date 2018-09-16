import numpy as np
import matplotlib.pyplot as plt
import json, glob
import sys, time, os
from numpy.linalg import inv
from unwrap.unwrap import unwrap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from linecache import getline, clearcache
from scipy.integrate import simps
from scipy.constants import *

from OCC.Display.SimpleGui import init_display
from OCC.gp import gp_Pnt, gp_Ax1, gp_Ax3, gp_Vec, gp_Dir
from OCC.gp import gp_Trsf, gp_Quaternion, gp_Pln
from OCC.TColgp  import TColgp_Array2OfPnt
from OCC.Geom    import Geom_BSplineSurface, Handle_Geom_BSplineSurface
from OCC.Geom    import Geom_Surface, Handle_Geom_Surface
from OCC.GeomAPI import GeomAPI_PointsToBSplineSurface
from OCC.GeomAPI import GeomAPI_ProjectPointOnSurf
from OCC.GeomAbs import GeomAbs_C2, GeomAbs_C0, GeomAbs_G1, GeomAbs_G2
from OCC.TopoDS  import TopoDS_Shape, TopoDS_Shell
from OCC.TopLoc  import TopLoc_Location
from OCC.BRep    import BRep_Tool_Surface, BRep_Builder
from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeFace
from OCC.BRepExtrema    import BRepExtrema_DistShapeShape
from OCCUtils.Construct import make_plane, make_line, make_wire, make_edge
from OCCUtils.Construct import vec_to_dir
from OCCUtils.Topology  import Topo, dumpTopology

def surf_trf (pnt, vec, face):
    axs = gp_Ax3 (pnt, vec_to_dir(vec))
    trf = gp_Trsf()
    trf.SetTransformation (axs, gp_Ax3())
    srf = face.Moved (TopLoc_Location (trf))
    return srf

def surf_spl(px, py, pz):
    nx, ny = px.shape
    pnt_2d = TColgp_Array2OfPnt(1, nx, 1, ny)
    for row in range (pnt_2d.LowerRow(), pnt_2d.UpperRow()+1):
        for col in range (pnt_2d.LowerCol(), pnt_2d.UpperCol()+1):
            i,j = row-1, col-1
            pnt = gp_Pnt (px[i, j], py[i, j], pz[i, j])
            pnt_2d.SetValue(row, col, pnt)
    
    curv = GeomAPI_PointsToBSplineSurface(pnt_2d, 3, 8, GeomAbs_G2, 0.001).Surface()
    surf = BRepBuilderAPI_MakeFace(curv, 1e-6).Face()
    return surf

def normal (geom, d0=0, d1=0):
    p, v0, v1 = gp_Pnt(), gp_Vec(), gp_Vec()
    geom.D1 (d0, d1, p, v0, v1)
    v0.Normalize()
    v1.Normalize()
    #v0.Reverse()
    #v1.Reverse()
    v2 = v1.Crossed (v0)
    v2.Normalize()
    print ("pnt", p.X() , p.Y() , p.Z() )
    print ("v_x", v0.X(), v0.Y(), v0.Z())
    print ("v_y", v1.X(), v1.Y(), v1.Z())
    print ("v_z", v2.X(), v2.Y(), v2.Z())
    return p, v0, v1, v2

def surf_dev (face, dx=50, dy=50):
    geom = BRep_Tool_Surface (face).GetObject()
    p0, v0_x, v0_y, v0_z = normal (geom, 0.0, 0.0)
    p1, v1_x, v1_y, v1_z = normal (geom, dx , 0.0)
    p2, v2_x, v2_y, v2_z = normal (geom, 0.0, dy )