import numpy as np
import matplotlib.pyplot as plt
import json
import glob
import sys
import time
import os
from linecache import getline, clearcache

from OCC.Display.SimpleGui import init_display
from OCC.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.gp import gp_Pln, gp_Trsf, gp_Lin
from OCC.Geom import Geom_Plane, Geom_Surface, Geom_BSplineSurface
from OCC.Geom import Geom_Curve, Geom_Line
from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeFace
from OCC.TColgp import TColgp_Array1OfPnt, TColgp_Array2OfPnt
from OCC.GeomAPI import GeomAPI_PointsToBSplineSurface, GeomAPI_IntCS
from OCC.GeomAPI import GeomAPI_ProjectPointOnSurf, GeomAPI_ProjectPointOnCurve
from OCC.GeomAbs import GeomAbs_C2, GeomAbs_C0, GeomAbs_G1, GeomAbs_G2
from OCC.GeomLProp import GeomLProp_SurfaceTool
from OCC.TopLoc import TopLoc_Location
from OCC.BRep import BRep_Tool
from OCCUtils.Construct import dir_to_vec, vec_to_dir
from OCCUtils.Construct import make_wire, make_edge, make_plane, make_line
from OCCUtils.Construct import project_edge_onto_plane, project_point_on_curve
from OCCUtils.Topology import Topo

from .ray_setup import get_axs, load_surface, reflect, axs_pln, get_deg
from ..pyocc.load import read_step_file
from ..pyocc.surface import surf_spl
from ..fileout import occ_to_grasp_cor, occ_to_grasp_rim


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


def make_edges(pts):
    edg = []
    for i in range(len(pts)-1):
        i0, i1 = i, i+1
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
    if os.path.exists(name+"_org.cor"):
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
        if os.path.exists(self.dir + self.name+"_org.cor"):
            self.axs_org = get_axs(self.dir + self.name + "_org.cor")
        else:
            self.axs_org = gp_Ax3()
        self.axs = get_axs(self.dir + self.name+".cor", self.local)
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


class RaySystem (object):

    def __init__(self, dir_name, ini_name, tar_name):
        self.dir = dir_name
        self.axs = gp_Ax3()
        self.ini = SurfSystem(self.dir, ini_name)
        self.tar = SurfSystem(self.dir, tar_name)
        self.display, self.start_display, self.add_menu, self.add_function_to_menu = init_display()

    def sxy_ini(self):
        prf = self.dir + self.ini.name + "_" + "sz" + "_profile.txt"
        sx = float(getline(prf, 7).split()[1])
        sy = float(getline(prf, 8).split()[1])
        sxy = gp_Pnt(sx, sy, 0).Transformed(self.ini.trf)
        h_srf = BRep_Tool.Surface(self.ini.srf)
        self.ini.sxy = GeomAPI_ProjectPointOnSurf(sxy, h_srf).Point(1)

    def sxy_tar(self):
        prf = self.dir + self.tar.name + "_" + "sz" + "_profile.txt"
        sx = float(getline(prf, 7).split()[1])
        sy = float(getline(prf, 8).split()[1])
        sxy = gp_Pnt(sx, sy, 0).Transformed(self.tar.trf)
        h_srf = BRep_Tool.Surface(self.tar.srf)
        self.tar.sxy = GeomAPI_ProjectPointOnSurf(sxy, h_srf).Point(1)

    def Reflect(self):
        h_surf = BRep_Tool.Surface(self.tar.srf)
        g_line = gp_Lin(self.ini.beam.Location(), self.ini.beam.Direction())
        ray = Geom_Line(g_line)
        if GeomAPI_IntCS(ray.GetHandle(), h_surf).NbPoints():
            self.tar.beam = reflect(self.ini.beam, self.tar.srf)
        else:
            pln = make_plane(
                self.tar.axs.Location(),
                dir_to_vec(self.tar.axs.Direction()),
                500, -500,
                500, -500
            )
            self.tar.beam = reflect(self.ini.beam, pln)

    def OptAxs(self, axs="y"):
        print(self.ini.name, "->", self.tar.name)
        if axs == None:
            ref = gp_Dir(0, 0, 1)
            pnt = gp_Pnt(0, 0, 0)
            ax = gp_Ax3(pnt, ref)
        elif axs == "x":
            ref = self.ini.beam.XDirection()
            pnt = self.tar.BeamLocal().Location()
            pnt.SetZ(0.0)
            pnt.SetY(0.0)
            ax = gp_Ax3(pnt, gp_Dir(0, 0, 1))
        elif axs == "y":
            ref = self.ini.beam.YDirection()
            pnt = self.tar.BeamLocal().Location()
            pnt.SetZ(0.0)
            pnt.SetX(0.0)
            ax = gp_Ax3(pnt, gp_Dir(0, 0, 1))
        else:
            ref = gp_Dir(0, 0, 1)
            pnt = gp_Pnt(0, 0, 0)
            ax = gp_Ax3(pnt, gp_Dir(0, 0, 1))
        print("Detect Target Position")
        print(self.ini.axs.Location())
        print(self.tar.axs.Location())
        print(self.ini.beam.Location())
        print(self.tar.beam.Location())
        print(self.ini.BeamLocal().Location())
        print(self.tar.BeamLocal().Location())

        ax = self.tar.Move_Axs(ax, gp_Ax3(), self.tar.axs)
        v0 = dir_to_vec(self.ini.beam.Direction())
        v1 = gp_Vec(self.ini.beam.Location(), ax.Location()).Normalized()
        rf = dir_to_vec(ref)
        deg = np.rad2deg(v0.AngleWithRef(v1, rf))

        print("Detect Rotation Angle")
        print(ax.Location())

        print(self.ini.name, axs, deg)
        print(v0)
        print(v1)
        print(rf)
        self.ini.RotAxs(deg/2, axs=axs)

    def Display_Shape(self, colors=["BLUE", "YELLOW"]):
        self.display.DisplayShape(axs_pln(gp_Ax3()))
        self.display.DisplayShape(axs_pln(self.ini.axs), color=colors[0])
        self.display.DisplayShape(axs_pln(self.tar.axs), color=colors[1])
        self.display.DisplayShape(self.ini.srf, color=colors[0])
        self.display.DisplayShape(self.tar.srf)
        self.display.DisplayShape(self.ini.beam.Location(), color=colors[0])
        self.display.DisplayShape(self.tar.beam.Location(), color=colors[1])
        self.display.DisplayShape(
            make_line(self.ini.beam.Location(), self.tar.beam.Location()), color=colors[0])
        self.display.FitAll()

    def Display_Update(self):
        self.display.EraseAll()
        self.display.DisplayShape(axs_pln(gp_Ax3()))
        self.display.DisplayShape(axs_pln(self.ini.axs))
        self.display.DisplayShape(axs_pln(self.tar.axs))
        self.display.DisplayShape(self.ini.srf)
        self.display.DisplayShape(self.tar.srf)
        self.display.DisplayShape(self.ini.sxy, color="BLUE")
        self.display.DisplayShape(self.tar.sxy, color="YELLOW")
        self.display.DisplayShape(make_line(self.ini.sxy, self.tar.sxy))
        self.display.FitAll()
        self.start_display()


class OptSystem (object):

    def __init__(self, dir_name, name, ini_name, tar_name):
        self.dir = dir_name
        self.axs = gp_Ax3()
        self.srf = SurfSystem(self.dir, name)
        self.ini = SurfSystem(self.dir, ini_name)
        self.tar = SurfSystem(self.dir, tar_name)
        self.display, self.start_display, self.add_menu, self.add_function_to_menu = init_display()

    def Reflect_ini(self):
        h_surf = BRep_Tool.Surface(self.srf.srf)
        g_line = gp_Lin(self.ini.beam.Location(), self.ini.beam.Direction())
        ray = Geom_Line(g_line)
        if GeomAPI_IntCS(ray.GetHandle(), h_surf).NbPoints():
            self.srf.beam = reflect(self.ini.beam, self.srf.srf)
        else:
            pln = make_plane(
                self.srf.axs.Location(),
                dir_to_vec(self.srf.axs.Direction()),
                500, -500,
                500, -500
            )
            self.srf.beam = reflect(self.ini.beam, pln)

    def Reflect_srf(self):
        h_surf = BRep_Tool.Surface(self.tar.srf)
        g_line = gp_Lin(self.srf.beam.Location(), self.srf.beam.Direction())
        ray = Geom_Line(g_line)
        if GeomAPI_IntCS(ray.GetHandle(), h_surf).NbPoints():
            self.tar.beam = reflect(self.srf.beam, self.tar.srf)
        else:
            pln = make_plane(
                self.tar.axs.Location(),
                dir_to_vec(self.tar.axs.Direction()),
                500, -500,
                500, -500
            )
            self.tar.beam = reflect(self.srf.beam, pln)

    def OptAxs(self, deg=0.1, axs="y"):
        print(self.ini.name, "->", self.srf.name, "->", self.tar.name)
        self.Reflect_ini()
        self.Reflect_srf()
        self.Display_Shape([None, "RED"])
        print(self.srf.beam.Location())
        print(self.tar.beam.Location())
        print(self.tar.axs.Location())
        ax = self.tar.GetTarget(axs)
        d0 = ax.Location().Distance(self.tar.BeamLocal().Location())
        d1 = d0

        for i in range(100):
            self.srf.RotAxs(deg, axs)
            self.Reflect_ini()
            self.Reflect_srf()

            ax = self.tar.GetTarget(axs)
            d0 = d1
            d1 = ax.Location().Distance(self.tar.BeamLocal().Location())
            print("Rot ", i, deg, d0, d1)

            if d1 > d0:
                deg = -deg/2
            elif d1 < 1/100:
                break
        print(self.srf.beam.Location())
        print(self.tar.beam.Location())
        print(self.tar.axs.Location())
        self.Display_Shape([None, "BLUE"])

    def Display_Shape(self, colors=["BLUE", "YELLOW"]):
        self.display.DisplayShape(axs_pln(gp_Ax3()))
        self.display.DisplayShape(axs_pln(self.srf.axs), color=colors[0])
        self.display.DisplayShape(axs_pln(self.ini.axs), color=colors[0])
        self.display.DisplayShape(axs_pln(self.tar.axs), color=colors[1])
        self.display.DisplayShape(self.ini.srf, color=colors[0])
        self.display.DisplayShape(self.srf.srf, color=colors[0])
        self.display.DisplayShape(self.tar.srf)
        self.display.DisplayShape(self.ini.beam.Location(), color=colors[0])
        self.display.DisplayShape(self.srf.beam.Location(), color=colors[0])
        self.display.DisplayShape(self.tar.beam.Location(), color=colors[1])
        self.display.DisplayShape(
            make_line(self.ini.beam.Location(), self.srf.beam.Location()), color=colors[0])
        self.display.DisplayShape(
            make_line(self.srf.beam.Location(), self.tar.beam.Location()), color=colors[0])
        self.display.FitAll()

    def Display_Update(self):
        self.display.EraseAll()
        self.display.DisplayShape(axs_pln(gp_Ax3()))
        self.display.DisplayShape(axs_pln(self.ini.axs))
        self.display.DisplayShape(axs_pln(self.tar.axs))
        self.display.DisplayShape(self.ini.srf)
        self.display.DisplayShape(self.tar.srf)
        self.display.DisplayShape(self.ini.sxy, color="BLUE")
        self.display.DisplayShape(self.tar.sxy, color="YELLOW")
        self.display.DisplayShape(make_line(self.ini.sxy, self.tar.sxy))
        self.display.FitAll()
        self.start_display()


class Multi_RaySystem (object):

    def __init__(self, dir_name, ini_name, tar_name):
        self.dir = dir_name
        self.axs = gp_Ax3()
        self.ini = SurfSystem(self.dir, ini_name)
        self.tar = SurfSystem(self.dir, tar_name)
        self.ini.MultiRay()
        self.tar.MultiRay()
        self.display, self.start_display, self.add_menu, self.add_function_to_menu = init_display()

    def Reflect(self):
        h_surf = BRep_Tool.Surface(self.tar.srf)
        g_line = gp_Lin(self.ini.beam.Location(), self.ini.beam.Direction())
        ray = Geom_Line(g_line)
        if GeomAPI_IntCS(ray.GetHandle(), h_surf).NbPoints():
            self.tar.beam = reflect(self.ini.beam, self.tar.srf)
        else:
            pln = make_plane(
                self.tar.axs.Location(),
                dir_to_vec(self.tar.axs.Direction()),
                500, -500,
                500, -500
            )
            self.tar.beam = reflect(self.ini.beam, pln)

    def BeamReflect(self, beam=gp_Ax3()):
        h_surf = BRep_Tool.Surface(self.tar.srf)
        g_line = gp_Lin(self.ini.beam.Location(), beam.Direction())
        ray = Geom_Line(g_line)
        if GeomAPI_IntCS(ray.GetHandle(), h_surf).NbPoints():
            tar_beam = reflect(beam, self.tar.srf)
        else:
            pln = make_plane(
                self.tar.axs.Location(),
                dir_to_vec(self.tar.axs.Direction()),
                500, -500,
                500, -500
            )
            tar_beam = reflect(beam, pln)
        return tar_beam

    def MultiReflect(self):
        self.tar.beam = self.BeamReflect(self.ini.beam)
        self.tar.beam_rght = self.BeamReflect(self.ini.beam_rght)
        self.tar.beam_left = self.BeamReflect(self.ini.beam_left)
        self.tar.beam_uppr = self.BeamReflect(self.ini.beam_uppr)
        self.tar.beam_bott = self.BeamReflect(self.ini.beam_bott)

    def OptAxs(self, axs="y"):
        print(self.ini.name, "->", self.tar.name)
        if axs == None:
            ref = gp_Dir(0, 0, 1)
            pnt = gp_Pnt(0, 0, 0)
            ax = gp_Ax3(pnt, ref)
        elif axs == "x":
            ref = self.ini.beam.XDirection()
            pnt = self.tar.BeamLocal().Location()
            pnt.SetZ(0.0)
            pnt.SetY(0.0)
            ax = gp_Ax3(pnt, gp_Dir(0, 0, 1))
        elif axs == "y":
            ref = self.ini.beam.YDirection()
            pnt = self.tar.BeamLocal().Location()
            pnt.SetZ(0.0)
            pnt.SetX(0.0)
            ax = gp_Ax3(pnt, gp_Dir(0, 0, 1))
        else:
            ref = gp_Dir(0, 0, 1)
            pnt = gp_Pnt(0, 0, 0)
            ax = gp_Ax3(pnt, gp_Dir(0, 0, 1))
        print("Detect Target Position")
        print(self.ini.axs.Location())
        print(self.tar.axs.Location())
        print(self.ini.beam.Location())
        print(self.tar.beam.Location())
        print(self.ini.BeamLocal().Location())
        print(self.tar.BeamLocal().Location())

        ax = self.tar.Move_Axs(ax, gp_Ax3(), self.tar.axs)
        v0 = dir_to_vec(self.ini.beam.Direction())
        v1 = gp_Vec(self.ini.beam.Location(), ax.Location()).Normalized()
        rf = dir_to_vec(ref)
        deg = np.rad2deg(v0.AngleWithRef(v1, rf))

        print("Detect Rotation Angle")
        print(ax.Location())

        print(self.ini.name, axs, deg)
        print(v0)
        print(v1)
        print(rf)
        self.ini.RotAxs(deg/2, axs=axs)

    def Display_Shape(self, colors=["BLUE", "YELLOW"]):
        self.display.DisplayShape(axs_pln(gp_Ax3()))
        self.display.DisplayShape(axs_pln(self.ini.axs), color=colors[0])
        self.display.DisplayShape(axs_pln(self.tar.axs), color=colors[1])
        self.display.DisplayShape(self.ini.srf, color=colors[0])
        self.display.DisplayShape(self.tar.srf)
        self.display.DisplayShape(self.ini.beam.Location(), color=colors[0])
        self.display.DisplayShape(self.tar.beam.Location(), color=colors[1])
        self.display.DisplayShape(
            make_line(self.ini.beam.Location(), self.tar.beam.Location()), color=colors[0])
        self.display.DisplayShape(
            make_line(self.ini.beam_rght.Location(), self.tar.beam_rght.Location()), color=colors[0])
        self.display.DisplayShape(
            make_line(self.ini.beam_left.Location(), self.tar.beam_left.Location()), color=colors[0])
        self.display.DisplayShape(
            make_line(self.ini.beam_uppr.Location(), self.tar.beam_uppr.Location()), color=colors[0])
        self.display.DisplayShape(
            make_line(self.ini.beam_bott.Location(), self.tar.beam_bott.Location()), color=colors[0])
        self.display.FitAll()

    def Display_Update(self):
        self.display.EraseAll()
        self.display.DisplayShape(axs_pln(gp_Ax3()))
        self.display.DisplayShape(axs_pln(self.ini.axs))
        self.display.DisplayShape(axs_pln(self.tar.axs))
        self.display.DisplayShape(self.ini.srf)
        self.display.DisplayShape(self.tar.srf)
        self.display.DisplayShape(self.ini.sxy, color="BLUE")
        self.display.DisplayShape(self.tar.sxy, color="YELLOW")
        self.display.DisplayShape(make_line(self.ini.sxy, self.tar.sxy))
        self.display.FitAll()
        self.start_display()


class Beam_RaySystem (object):

    def __init__(self, dir_name, ini_name, tar_name):
        self.dir = dir_name
        self.axs = gp_Ax3()
        self.ini = SurfSystem(self.dir, ini_name)
        self.tar = SurfSystem(self.dir, tar_name)
        self.ini.MultiRay()
        self.tar.MultiRay()
        self.ini.ph_0 = wavefront(rxy=[0, 0], axs=self.ini.beam)
        self.display, self.start_display, self.add_menu, self.add_function_to_menu = init_display()

    def Reflect(self):
        h_surf = BRep_Tool.Surface(self.tar.srf)
        g_line = gp_Lin(self.ini.beam.Location(), self.ini.beam.Direction())
        ray = Geom_Line(g_line)
        if GeomAPI_IntCS(ray.GetHandle(), h_surf).NbPoints():
            self.tar.beam = reflect(self.ini.beam, self.tar.srf)
        else:
            pln = make_plane(
                self.tar.axs.Location(),
                dir_to_vec(self.tar.axs.Direction()),
                500, -500,
                500, -500
            )
            self.tar.beam = reflect(self.ini.beam, pln)

    def BeamReflect(self, beam=gp_Ax3()):
        h_surf = BRep_Tool.Surface(self.tar.srf)
        g_line = gp_Lin(self.ini.beam.Location(), beam.Direction())
        ray = Geom_Line(g_line)
        if GeomAPI_IntCS(ray.GetHandle(), h_surf).NbPoints():
            tar_beam = reflect(beam, self.tar.srf)
        else:
            pln = make_plane(
                self.tar.axs.Location(),
                dir_to_vec(self.tar.axs.Direction()),
                500, -500,
                500, -500
            )
            tar_beam = reflect(beam, pln)
        return tar_beam

    def MultiReflect(self):
        self.tar.beam = self.BeamReflect(self.ini.beam)
        self.tar.beam_rght = self.BeamReflect(self.ini.beam_rght)
        self.tar.beam_left = self.BeamReflect(self.ini.beam_left)
        self.tar.beam_uppr = self.BeamReflect(self.ini.beam_uppr)
        self.tar.beam_bott = self.BeamReflect(self.ini.beam_bott)

    def OptAxs(self, axs="y"):
        print(self.ini.name, "->", self.tar.name)
        if axs == None:
            ref = gp_Dir(0, 0, 1)
            pnt = gp_Pnt(0, 0, 0)
            ax = gp_Ax3(pnt, ref)
        elif axs == "x":
            ref = self.ini.beam.XDirection()
            pnt = self.tar.BeamLocal().Location()
            pnt.SetZ(0.0)
            pnt.SetY(0.0)
            ax = gp_Ax3(pnt, gp_Dir(0, 0, 1))
        elif axs == "y":
            ref = self.ini.beam.YDirection()
            pnt = self.tar.BeamLocal().Location()
            pnt.SetZ(0.0)
            pnt.SetX(0.0)
            ax = gp_Ax3(pnt, gp_Dir(0, 0, 1))
        else:
            ref = gp_Dir(0, 0, 1)
            pnt = gp_Pnt(0, 0, 0)
            ax = gp_Ax3(pnt, gp_Dir(0, 0, 1))
        print("Detect Target Position")
        print(self.ini.axs.Location())
        print(self.tar.axs.Location())
        print(self.ini.beam.Location())
        print(self.tar.beam.Location())
        print(self.ini.BeamLocal().Location())
        print(self.tar.BeamLocal().Location())

        ax = self.tar.Move_Axs(ax, gp_Ax3(), self.tar.axs)
        v0 = dir_to_vec(self.ini.beam.Direction())
        v1 = gp_Vec(self.ini.beam.Location(), ax.Location()).Normalized()
        rf = dir_to_vec(ref)
        deg = np.rad2deg(v0.AngleWithRef(v1, rf))

        print("Detect Rotation Angle")
        print(ax.Location())

        print(self.ini.name, axs, deg)
        print(v0)
        print(v1)
        print(rf)
        self.ini.RotAxs(deg/2, axs=axs)

    def Display_Shape(self, colors=["BLUE", "YELLOW"]):
        self.display.DisplayShape(axs_pln(gp_Ax3()))
        self.display.DisplayShape(axs_pln(self.ini.axs), color=colors[0])
        self.display.DisplayShape(axs_pln(self.tar.axs), color=colors[1])
        self.display.DisplayShape(self.ini.srf, color=colors[0])
        self.display.DisplayShape(self.tar.srf)
        self.display.DisplayShape(self.ini.beam.Location(), color=colors[0])
        self.display.DisplayShape(self.tar.beam.Location(), color=colors[1])
        self.display.DisplayShape(
            make_line(self.ini.beam.Location(), self.tar.beam.Location()), color=colors[0])
        self.display.DisplayShape(
            make_line(self.ini.beam_rght.Location(), self.tar.beam_rght.Location()), color=colors[0])
        self.display.DisplayShape(
            make_line(self.ini.beam_left.Location(), self.tar.beam_left.Location()), color=colors[0])
        self.display.DisplayShape(
            make_line(self.ini.beam_uppr.Location(), self.tar.beam_uppr.Location()), color=colors[0])
        self.display.DisplayShape(
            make_line(self.ini.beam_bott.Location(), self.tar.beam_bott.Location()), color=colors[0])
        self.display.FitAll()

    def Display_Update(self):
        self.display.EraseAll()
        self.display.DisplayShape(axs_pln(gp_Ax3()))
        self.display.DisplayShape(axs_pln(self.ini.axs))
        self.display.DisplayShape(axs_pln(self.tar.axs))
        self.display.DisplayShape(self.ini.srf)
        self.display.DisplayShape(self.tar.srf)
        self.display.DisplayShape(self.ini.sxy, color="BLUE")
        self.display.DisplayShape(self.tar.sxy, color="YELLOW")
        self.display.DisplayShape(make_line(self.ini.sxy, self.tar.sxy))
        self.display.FitAll()
        self.start_display()
