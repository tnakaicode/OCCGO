import numpy as np
import matplotlib.pyplot as plt
import json
import glob
import sys
import time
import os
from linecache import getline, clearcache

from OCC.Display.SimpleGui import init_display
from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_Pln, gp_Trsf, gp_Lin
from OCC.Core.Geom import Geom_Plane, Geom_Surface, Geom_BSplineSurface
from OCC.Core.Geom import Geom_Curve, Geom_Line
from OCC.Core.GeomLProp import GeomLProp_SurfaceTool, GeomLProp_SLProps
from OCC.Core.GeomAbs import GeomAbs_C2, GeomAbs_C0, GeomAbs_G1, GeomAbs_G2
from OCC.Core.GeomAPI import GeomAPI_ProjectPointOnSurf, GeomAPI_ProjectPointOnCurve
from OCC.Core.GeomAPI import GeomAPI_PointsToBSplineSurface, GeomAPI_IntCS
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeFace
from OCC.Core.TColgp import TColgp_Array1OfPnt, TColgp_Array2OfPnt
from OCC.Core.TopLoc import TopLoc_Location
from OCC.Core.BRep import BRep_Tool
from OCCUtils.Construct import dir_to_vec, vec_to_dir
from OCCUtils.Construct import make_wire, make_edge, make_plane, make_line
from OCCUtils.Construct import project_edge_onto_plane, project_point_on_curve
from OCCUtils.Topology import Topo

from .ray_setup import get_axs, load_surface, reflect, axs_pln, get_deg
from ..pyocc.load import read_step_file
from ..pyocc.surface import surf_spl
from ..fileout import occ_to_grasp_cor, occ_to_grasp_rim
from ..geomtory import curvature, fit_surf


def wavefront(rxy=[1000, 1000], axs=gp_Ax3()):
    px = np.linspace(-1, 1, 100) * 50
    py = np.linspace(-1, 1, 100) * 50
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


def axs_curvature(h_surf, u=0, v=0):
    prop = GeomLProp_SLProps(2, 0.01)
    prop.SetSurface(h_surf)
    prop.SetParameters(u, v)

    d1, d2 = gp_Dir(), gp_Dir()
    prop.CurvatureDirections(d1, d2)
    vz = dir_to_vec(prop.Normal())
    v1 = dir_to_vec(d1)
    v2 = dir_to_vec(d2)
    c1 = prop.MaxCurvature()
    c2 = prop.MinCurvature()

    if c1 == 0:
        r1 = 0
    else:
        r1 = 1 / c1

    if c2 == 0:
        r2 = 0
    else:
        r2 = 1 / c2

    print("Max", c1, r1, v1)
    print("Min", c2, r1, v2)
    print(v1.Dot(v2))
    print(prop.Value())
    return vz, v1, v2, r1, r2


def make_edges(pts):
    edg = []
    for i in range(len(pts) - 1):
        i0, i1 = i, i + 1
        edg.append(make_edge(pts[i0], pts[i1]))
    return make_wire(edg)


def set_surface(filename):
    if os.path.exists(filename):
        shp = read_step_file(filename)
        for face in Topo(shp).faces():
            surf = face
    else:
        surf = make_plane()
    return surf


def set_axs_pln(name):
    if os.path.exists(name + "_org.cor"):
        ax = get_axs(name + "_org.cor")
    else:
        ax = gp_Ax3()
    axs = get_axs(name + ".cor", ax)
    srf = set_surface(name + ".stp")
    trf = gp_Trsf()
    trf.SetTransformation(axs, gp_Ax3())
    loc_face = TopLoc_Location(trf)
    srf.Location(loc_face)
    return axs, srf, trf


def setup_sxy(filename, ext="ex"):
    axs, srf, trf = set_axs_pln(filename)
    prf = filename + "_" + ext + "_profile.txt"
    sx = float(getline(prf, 7).split()[1])
    sy = float(getline(prf, 8).split()[1])
    sxy = gp_Pnt(sx, sy, 0).Transformed(trf)
    h_srf = BRep_Tool.Surface(srf)
    sxy = GeomAPI_ProjectPointOnSurf(sxy, h_srf).Point(1)
    return axs, srf, trf, sxy


class SurfSystem (object):

    def __init__(self, dir_name, name, local=gp_Ax3()):
        self.dir = dir_name
        self.name = name
        self.local = local
        self.Import()
        self.beam = self.axs

    def GetTarget(self, axs="y"):
        if axs == None:
            pnt = gp_Pnt(0, 0, 0)
            ax = gp_Ax3(pnt, gp_Dir(0, 0, 1), self.axs.XDirection())
        elif axs == "x":
            pnt = self.BeamLocal().Location()
            pnt.SetZ(0.0)
            pnt.SetY(0.0)
            ax = gp_Ax3(pnt, gp_Dir(0, 0, 1), self.axs.XDirection())
        elif axs == "y":
            pnt = self.BeamLocal().Location()
            pnt.SetZ(0.0)
            pnt.SetX(0.0)
            ax = gp_Ax3(pnt, gp_Dir(0, 0, 1), self.axs.XDirection())
        else:
            ref = gp_Dir(0, 0, 1)
            pnt = gp_Pnt(0, 0, 0)
            ax = gp_Ax3(pnt, gp_Dir(0, 0, 1), self.axs.XDirection())
        return ax

    def Import(self):
        if os.path.exists(self.dir + self.name + "_org.cor"):
            self.axs_org = get_axs(self.dir + self.name + "_org.cor")
        else:
            self.axs_org = gp_Ax3()
        self.axs = get_axs(self.dir + self.name + ".cor", self.local)
        self.srf = set_surface(self.dir + self.name + ".stp")
        self.trf = gp_Trsf()
        self.trf.SetTransformation(self.axs, self.local)
        self.loc_face = TopLoc_Location(self.trf)
        self.srf.Location(self.loc_face)

    def Move_Beam(self, ax0=gp_Ax3(), ax1=gp_Ax3()):
        trf = gp_Trsf()
        trf.SetTransformation(ax1, ax0)
        loc = TopLoc_Location(trf)
        self.beam.Transform(trf)

    def Move_Axs(self, axs=gp_Ax3(), ax0=gp_Ax3(), ax1=gp_Ax3()):
        trf = gp_Trsf()
        trf.SetTransformation(ax1, ax0)
        loc = TopLoc_Location(trf)
        return axs.Transformed(trf)

    def RotAxs(self, deg, axs="y"):
        if axs == None:
            rot_axs = gp_Ax1(self.axs.Location(), self.axs.XDirection())
        elif axs == "x":
            rot_axs = gp_Ax1(self.axs.Location(), self.axs.XDirection())
        elif axs == "y":
            rot_axs = gp_Ax1(self.axs.Location(), self.axs.YDirection())
        else:
            rot_axs = gp_Ax1(self.axs.Location(), self.axs.YDirection())
        self.Move_Beam(gp_Ax3(), gp_Ax3().Rotated(rot_axs, np.deg2rad(deg)))
        self.axs.Rotate(rot_axs, np.deg2rad(deg))
        self.srf = set_surface(self.dir + self.name + ".stp")
        self.trf = gp_Trsf()
        self.trf.SetTransformation(self.axs, gp_Ax3())
        self.loc_face = TopLoc_Location(self.trf)
        self.srf.Location(self.loc_face)

    def BeamLocal(self):
        return self.Move_Axs(self.beam, self.axs, gp_Ax3())

    def BeamSave(self):
        get_deg(self.axs, dir_to_vec(self.beam.Direction()))
        occ_to_grasp_cor(self.BeamLocal(), name=self.name,
                         filename=self.dir + self.name + "_beam.cor")

    def AxsLocal(self):
        return self.Move_Axs(self.axs, self.axs, self.local)

    def MultiRay(self, rght=[-10, 0], left=[10, 0], uppr=[0, 10], bott=[0, -10]):
        ax_rght = gp_Ax3(gp_Pnt(*rght, 0), gp_Dir(0, 0, 1))
        ax_left = gp_Ax3(gp_Pnt(*left, 0), gp_Dir(0, 0, 1))
        ax_uppr = gp_Ax3(gp_Pnt(*uppr, 0), gp_Dir(0, 0, 1))
        ax_bott = gp_Ax3(gp_Pnt(*bott, 0), gp_Dir(0, 0, 1))
        self.beam_rght = self.Move_Axs(self.beam, gp_Ax3(), ax_rght)
        self.beam_left = self.Move_Axs(self.beam, gp_Ax3(), ax_left)
        self.beam_uppr = self.Move_Axs(self.beam, gp_Ax3(), ax_uppr)
        self.beam_bott = self.Move_Axs(self.beam, gp_Ax3(), ax_bott)


class GaussSystem (SurfSystem):

    def __init__(self, dir_name, name, local=gp_Ax3()):
        super(GaussSystem, self).__init__(dir_name, name, local)

    def Init_Beam(self):
        prf = self.dir + self.name + "_ex_profile.txt"
        sx = float(getline(prf, 7).split()[1])
        sy = float(getline(prf, 8).split()[1])
        wx = float(getline(prf, 9).split()[1])
        wy = float(getline(prf, 10).split()[1])
        rt = float(getline(prf, 11).split()[1])
        dx = float(getline(prf, 12).split()[1])
        dy = float(getline(prf, 13).split()[1])
        rx = float(getline(prf, 14).split()[1])
        ry = float(getline(prf, 15).split()[1])
        print(prf)
        print(sx, sy)
        print(wx, wy)
        print(rt)
        print(dx, dy)
        print(rx, ry)

        self.wave = wavefront([-rx, -ry], self.beam)
