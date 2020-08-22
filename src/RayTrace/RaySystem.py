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
from .SurfSystem import SurfSystem, GaussSystem, wavefront, axs_curvature
from .SurfSystem import make_edges, set_surface, set_axs_pln, setup_sxy
from ..pyocc.load import read_step_file
from ..pyocc.surface import surf_spl
from ..fileout import occ_to_grasp_cor, occ_to_grasp_rim
from ..geomtory import curvature, fit_surf


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
        if GeomAPI_IntCS(ray, h_surf).NbPoints():
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
        self.ini.RotAxs(deg / 2, axs=axs)

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
        if GeomAPI_IntCS(ray, h_surf).NbPoints():
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
        if GeomAPI_IntCS(ray, h_surf).NbPoints():
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
                deg = -deg / 2
            elif d1 < 1 / 100:
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
        if GeomAPI_IntCS(ray, h_surf).NbPoints():
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
        if GeomAPI_IntCS(ray, h_surf).NbPoints():
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
        self.ini.RotAxs(deg / 2, axs=axs)

    def Display_Shape(self, colors=["BLUE", "YELLOW"]):
        #self.display.DisplayShape(axs_pln(gp_Ax3()))
        #self.display.DisplayShape(axs_pln(self.ini.axs), color=colors[0])
        #self.display.DisplayShape(axs_pln(self.tar.axs), color=colors[1])
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
        if GeomAPI_IntCS(ray, h_surf).NbPoints():
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
        if GeomAPI_IntCS(ray, h_surf).NbPoints():
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
        self.ini.RotAxs(deg / 2, axs=axs)

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


class GOSystem (object):

    def __init__(self, dir_name, ini_name, tar_name, wave=1.0):
        self.dir = dir_name
        self.axs = gp_Ax3()
        self.ini = GaussSystem(dir_name, ini_name)
        self.tar = GaussSystem(dir_name, tar_name)
        self.display, self.start_display, self.add_menu, self.add_function_to_menu = init_display()

        self.wave = wave
        self.knum = 2 * np.pi / wave

    def Reflect(self):
        h_surf = BRep_Tool.Surface(self.tar.srf)
        g_line = gp_Lin(self.ini.beam.Location(), self.ini.beam.Direction())
        ray = Geom_Line(g_line)
        if GeomAPI_IntCS(ray, h_surf).NbPoints():
            self.tar.beam = reflect(self.ini.beam, self.tar.srf)
        else:
            pln = make_plane(
                self.tar.axs.Location(),
                dir_to_vec(self.tar.axs.Direction()),
                500, -500,
                500, -500
            )
            h_surf = BRep_Tool.Surface(pln)
            self.tar.beam = reflect(self.ini.beam, pln)

        print(self.ini.beam.Location())
        print(self.tar.beam.Location())
        
        GeomAPI_IntCS(ray, h_surf).IsDone()
        uvw = GeomAPI_IntCS(ray, h_surf).Parameters(1)
        u, v, w = uvw
        print(u, v, w)
        vz, v1, v2, r1, r2 = axs_curvature(h_surf, u, v)

        tar_surf_axs = gp_Ax3(self.tar.beam.Location(),
                              vec_to_dir(vz), vec_to_dir(v1))
        tar_surf = wavefront([r1, r2], tar_surf_axs)
        self.display.DisplayShape(tar_surf, color="BLUE")
        self.display.DisplayShape(axs_pln(tar_surf_axs))

        self.GO_Prop(w)

        h_tar_wave = BRep_Tool.Surface(self.ini_wave)
        vz, v1, v2, r1, r2 = axs_curvature(h_tar_wave, u, v)
        tar_wave_axs = self.tar.beam.Translated(gp_Vec(0,0,0))
        tar_wave_axs.SetXDirection(vec_to_dir(v1))
        self.tar.wave = wavefront([r1, r2], tar_wave_axs)
        self.display.DisplayShape(self.tar.wave,  color="RED")
        self.display.DisplayShape(axs_pln(tar_wave_axs))

    def GO_Prop(self, s=0):
        h_ini_wave = BRep_Tool.Surface(self.ini.wave)
        vz, v1, v2, r1, r2 = axs_curvature(h_ini_wave, 0.5, 0.5)
        
        r1_z = r1 + s / 2
        r2_z = r2 + s / 2
        ini_wave_axs = self.ini.beam.Translated(
            dir_to_vec(self.ini.beam.Direction()).Scaled(s / 2))
        ini_wave_axs.SetXDirection(vec_to_dir(v1))
        ini_wave = wavefront([r1_z, r2_z], ini_wave_axs)
        self.display.DisplayShape(ini_wave)

        r1_z = r1 + s
        r2_z = r2 + s
        self.ini_wave_axs = self.ini.beam.Translated(
            dir_to_vec(self.ini.beam.Direction()).Scaled(s))
        self.ini_wave_axs.SetXDirection(vec_to_dir(v1))
        self.ini_wave = wavefront([r1_z, r2_z], self.ini_wave_axs)
        self.display.DisplayShape(self.ini_wave)

    def BeamReflect(self, beam=gp_Ax3()):
        h_surf = BRep_Tool.Surface(self.tar.srf)
        g_line = gp_Lin(self.ini.beam.Location(), beam.Direction())
        ray = Geom_Line(g_line)
        if GeomAPI_IntCS(ray, h_surf).NbPoints():
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
        self.ini.RotAxs(deg / 2, axs=axs)

    def Display_Shape(self, colors=["BLUE", "YELLOW"]):
        self.display.DisplayShape(axs_pln(gp_Ax3()))
        self.display.DisplayShape(axs_pln(self.ini.axs), color=colors[0])
        self.display.DisplayShape(axs_pln(self.tar.axs), color=colors[1])
        self.display.DisplayShape(self.ini.srf, color=colors[0])
        self.display.DisplayShape(self.tar.srf, color=colors[1])
        self.display.DisplayShape(
            make_line(self.ini.beam.Location(), self.tar.beam.Location()))
        self.display.DisplayShape(self.ini.beam.Location(), color=colors[0])
        self.display.DisplayShape(self.tar.beam.Location(), color=colors[1])
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
