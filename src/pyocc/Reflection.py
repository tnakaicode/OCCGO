import numpy as np
from linecache import getline, clearcache

from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_Pln, gp_Lin
from OCC.Core.BRep import BRep_Tool
from OCC.Core.Geom import Geom_Plane, Geom_Surface, Geom_BSplineSurface
from OCC.Core.Geom import Geom_Curve, Geom_Line, Geom_Ellipse
from OCC.Core.GeomAPI import GeomAPI_PointsToBSplineSurface, GeomAPI_IntCS
from OCC.Core.GeomAbs import GeomAbs_C2, GeomAbs_C0, GeomAbs_G1, GeomAbs_G2
from OCC.Core.GeomLProp import GeomLProp_SurfaceTool
from OCC.Core.TColgp import TColgp_Array1OfPnt, TColgp_Array2OfPnt
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeFace
from OCCUtils.Construct import dir_to_vec, vec_to_dir
from OCCUtils.Construct import project_edge_onto_plane
from OCC.Core.GeomProjLib import geomprojlib_Project


def reflect(p0, v0, face):
    h_surf = BRep_Tool.Surface(face)
    ray = Geom_Line(gp_Lin(p0, vec_to_dir(v0)))
    uvw = GeomAPI_IntCS(ray.GetHandle(), h_surf).Parameters(1)
    u, v, w = uvw
    p1, vx, vy = gp_Pnt(), gp_Vec(), gp_Vec()
    GeomLProp_SurfaceTool.D1(h_surf, u, v, p1, vx, vy)
    vz = vx.Crossed(vy)
    vx.Normalize()
    vy.Normalize()
    vz.Normalize()
    v1 = v0.Mirrored(gp_Ax2(p1, vec_to_dir(vz)))
    return p1, v1


def get_surface_sfc(filename):
    nx, ny = [int(s) for s in getline(filename, 3).split()]
    xs, ys, xe, ye = [float(s) for s in getline(filename, 2).split()]
    px = np.linspace(xs, xe, nx)
    py = np.linspace(ys, ye, ny)
    mesh = np.meshgrid(px, py)
    surf = np.loadtxt(filename, skiprows=3).T
    nx, ny = surf.shape
    pnt_2d = TColgp_Array2OfPnt(1, nx, 1, ny)
    for row in range(pnt_2d.LowerRow(), pnt_2d.UpperRow()+1):
        for col in range(pnt_2d.LowerCol(), pnt_2d.UpperCol()+1):
            i, j = row-1, col-1
            pnt = gp_Pnt(mesh[0][i, j], mesh[1][i, j], surf[i, j])
            pnt_2d.SetValue(row, col, pnt)
    surface = GeomAPI_PointsToBSplineSurface(
        pnt_2d, 3, 8, GeomAbs_G2, 0.001).Surface()
    srf_face = BRepBuilderAPI_MakeFace(surface, 0, 1, 0, 1, 0.001).Face()
    return srf_face


def get_ellips(axs, wxy, face=None):
    wx, wy = wxy
    if wx >= wy:
        ax2 = axs
        r0, r1 = wx, wy
    else:
        ax2 = axs.Rotated(axs.Axis(), np.deg2rad(90))
        r0, r1 = wy, wx
    el = Geom_Ellipse(ax2, r0, r1)

    if face == None:
        return el
    else:
        curv = geomprojlib_Project(el.GetHandle(), BRep_Tool.Surface(face))
        return curv


def Prj_pnt_to_face(axs, pnt, face):
    lin = gp_Lin(pnt, axs.Direction())
    sxy = GeomAPI_IntCS(Geom_Line(lin).GetHandle(),
                        BRep_Tool.Surface(face)).Point(1)
    return sxy
