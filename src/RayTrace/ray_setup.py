import numpy as np
import matplotlib.pyplot as plt
import json
import glob
import sys
import time
import os

from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir, gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_Pln, gp_Trsf, gp_Lin
from OCC.Core.Geom import Geom_Plane, Geom_Surface, Geom_BSplineSurface
from OCC.Core.Geom import Geom_Curve, Geom_Line
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeFace
from OCC.Core.TColgp import TColgp_Array1OfPnt, TColgp_Array2OfPnt
from OCC.Core.GeomAPI import GeomAPI_PointsToBSplineSurface, GeomAPI_IntCS
from OCC.Core.GeomAPI import GeomAPI_ProjectPointOnSurf, GeomAPI_ProjectPointOnCurve
from OCC.Core.GeomAbs import GeomAbs_C2, GeomAbs_C0, GeomAbs_G1, GeomAbs_G2
from OCC.Core.GeomLProp import GeomLProp_SurfaceTool, GeomLProp_SLProps
from OCC.Core.TopLoc import TopLoc_Location
from OCC.Core.BRep import BRep_Tool
from OCC.Core.TopoDS import TopoDS_Shape
from OCCUtils.Construct import dir_to_vec, vec_to_dir
from OCCUtils.Construct import make_wire, make_edge, make_plane, make_line
from OCCUtils.Construct import project_edge_onto_plane, project_point_on_curve
from OCCUtils.Topology import Topo


def load_surface(mesh, surf):
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


def reflect(beam, face, axs=gp_Ax3()):
    p0, v0 = beam.Location(), dir_to_vec(beam.Direction())
    h_surf = BRep_Tool.Surface(face)
    ray = Geom_Line(gp_Lin(p0, vec_to_dir(v0)))
    if GeomAPI_IntCS(ray.GetHandle(), h_surf).NbPoints() == 0:
        print("Out of Surface", axs.Location())
        pln = make_plane(
            axs.Location(), dir_to_vec(axs.Direction()), 500, -500, 500, -500
        )
        h_surf = BRep_Tool.Surface(pln)
    GeomAPI_IntCS(ray.GetHandle(), h_surf).IsDone()
    uvw = GeomAPI_IntCS(ray.GetHandle(), h_surf).Parameters(1)
    u, v, w = uvw
    p1, vx, vy = gp_Pnt(), gp_Vec(), gp_Vec()
    GeomLProp_SurfaceTool.D1(h_surf, u, v, p1, vx, vy)
    vz = vx.Crossed(vy)
    vx.Normalize()
    vy.Normalize()
    vz.Normalize()
    v1 = v0.Mirrored(gp_Ax2(p1, vec_to_dir(vz), vec_to_dir(vx)))
    return gp_Ax3(p1, vec_to_dir(v1), beam.XDirection().Reversed())


def axs_pln(axs):
    pnt = axs.Location()
    vx = dir_to_vec(axs.XDirection()).Scaled(10)
    vy = dir_to_vec(axs.YDirection()).Scaled(20)
    vz = dir_to_vec(axs.Direction()).Scaled(30)
    lx = make_line(pnt, gp_Pnt((gp_Vec(pnt.XYZ())+vx).XYZ()))
    ly = make_line(pnt, gp_Pnt((gp_Vec(pnt.XYZ())+vy).XYZ()))
    lz = make_line(pnt, gp_Pnt((gp_Vec(pnt.XYZ())+vz).XYZ()))
    return lx, ly, lz


def Transform_obj(axs, obj):
    trf = gp_Trsf()
    trf.SetTransformation(gp_Ax3(), axs)


def Transform_axs(ax0=gp_Ax3(), ax1=gp_Ax3()):
    trf = gp_Trsf()
    trf.SetTransformation(gp_Ax3(), ax1)
    ax0.Transform(trf)


def get_axs(filename, ax=gp_Ax3()):
    dat = np.loadtxt(filename, skiprows=2)
    pnt = gp_Pnt(*dat[0])
    d_x = gp_Dir(*dat[1])
    d_y = gp_Dir(*dat[2])
    d_z = d_x.Crossed(d_y)
    axs = gp_Ax3(pnt, d_z, d_x)
    trf = gp_Trsf()
    trf.SetTransformation(ax, gp_Ax3())
    axs.Transform(trf)
    return axs


def get_deg(axs, vec):
    vx = dir_to_vec(axs.XDirection())
    vy = dir_to_vec(axs.YDirection())
    vz = dir_to_vec(axs.Direction())
    pln_x = Geom_Plane(axs.Location(), axs.YDirection())
    pln_y = Geom_Plane(axs.Location(), axs.XDirection())
    vec_p = gp_Pnt((gp_Vec(axs.Location().XYZ()) + vec).XYZ())
    pnt_x = GeomAPI_ProjectPointOnSurf(vec_p, pln_x.GetHandle()).Point(1)
    pnt_y = GeomAPI_ProjectPointOnSurf(vec_p, pln_y.GetHandle()).Point(1)
    vec_x = gp_Vec(axs.Location(), pnt_x)
    vec_y = gp_Vec(axs.Location(), pnt_y)
    deg_x = vec_x.AngleWithRef(vz, vy)
    deg_y = vec_y.AngleWithRef(vz, vx)
    print(np.rad2deg(deg_x), np.rad2deg(deg_y))
    return deg_x, deg_y
