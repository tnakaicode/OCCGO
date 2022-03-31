import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
from linecache import getline, clearcache
import argparse

sys.path.append(os.path.join("../../"))
from src.base_occ import dispocc, set_trf
from src.geomtory import curvature, grasp_sfc

import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)

from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_Pln, gp_Circ
from OCC.Core.gp import gp_Trsf, gp_Quaternion
from OCC.Core.TopoDS import TopoDS_Shape, TopoDS_Compound
from OCC.Core.TopLoc import TopLoc_Location
from OCC.Core.TColgp import TColgp_Array2OfPnt
from OCC.Core.Geom import Geom_Curve, Geom_Plane, Geom_Line
from OCC.Core.GeomAbs import GeomAbs_C2, GeomAbs_C0, GeomAbs_G1, GeomAbs_G2
from OCC.Core.GeomAPI import GeomAPI_PointsToBSplineSurface
from OCC.Core.GeomAPI import GeomAPI_ProjectPointOnSurf
from OCC.Core.GeomAPI import GeomAPI_IntCS, GeomAPI_IntSS
from OCC.Core.GeomLProp import GeomLProp_SurfaceTool, GeomLProp_CurveTool
from OCC.Core.BRep import BRep_Tool
from OCC.Core.BRepTools import breptools
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeFace
from OCC.Core.BRepProj import BRepProj_Projection
from OCCUtils.Topology import Topo
from OCCUtils.Construct import make_box, make_line, make_wire, make_edge
from OCCUtils.Construct import make_plane, make_polygon
from OCCUtils.Construct import point_to_vector, vector_to_point
from OCCUtils.Construct import dir_to_vec, vec_to_dir


def float_to_string(number):
    if number == 0 or abs(np.log10(abs(number))) < 100:
        return ' {: 0.10E}'.format(number)
    else:
        return ' {: 0.10E}'.format(number).replace('E', '')


def surf_spl_pcd(px, py, pz):
    nx, ny = px.shape
    pnt_2d = TColgp_Array2OfPnt(1, nx, 1, ny)
    for row in range(pnt_2d.LowerRow(), pnt_2d.UpperRow() + 1):
        for col in range(pnt_2d.LowerCol(), pnt_2d.UpperCol() + 1):
            i, j = row - 1, col - 1
            pnt = gp_Pnt(px[i, j], py[i, j], pz[i, j])
            pnt_2d.SetValue(row, col, pnt)
    curv = GeomAPI_PointsToBSplineSurface(
        pnt_2d, 3, 8, GeomAbs_G2, 0.001).Surface()
    surf = BRepBuilderAPI_MakeFace(curv, 1e-6).Face()
    return surf, pnt_2d


class SurfOCC(object):

    def __init__(self, axs=gp_Ax3()):
        self.rot = axs
        self.axs = gp_Ax3(self.rot.Ax2())
        self.rim = make_edge(gp_Circ(self.axs.Ax2(), 100))
        self.pln = dispocc.make_plane_axs(self.axs)
        self.surf = make_plane(
            self.axs.Location(), dir_to_vec(self.axs.Direction()),
            -500, 500, -500, 500)

    def get_trsf(self):
        self.trf = gp_Trsf()
        self.trf.SetTransformation(self.axs, gp_Ax3())
        return self.trf



    def RotateSurf(self, deg=0.0, axs="z"):
        if axs == "x":
            ax1 = gp_Ax1(self.rot.Location(), self.rot.XDirection())
        elif axs == "y":
            ax1 = gp_Ax1(self.rot.Location(), self.rot.YDirection())
        elif axs == "z":
            ax1 = gp_Ax1(self.rot.Location(), self.rot.Direction())
        else:
            ax1 = gp_Ax1(self.rot.Location(), self.rot.Direction())
        trf = gp_Trsf()
        trf.SetRotation(ax1, np.deg2rad(deg))
        self.rot.Transform(trf)
        self.axs.Transform(trf)
        self.surf.Move(TopLoc_Location(trf))

    def MovXYZSurf(self, dst=0.0, axs="z"):
        if axs == "x":
            vec = dir_to_vec(self.rot.XDirection())
        elif axs == "y":
            vec = dir_to_vec(self.rot.YDirection())
        elif axs == "z":
            vec = dir_to_vec(self.rot.Direction())
        else:
            vec = dir_to_vec(self.rot.Direction())
        trf = gp_Trsf()
        trf.SetTranslation(vec.Scaled(dst))
        self.rot.Transform(trf)
        self.axs.Transform(trf)
        self.surf.Move(TopLoc_Location(trf))

    def RotateAxis(self, deg=0.0, axs="z"):
        if axs == "x":
            ax1 = gp_Ax1(self.rot.Location(), self.rot.XDirection())
        elif axs == "y":
            ax1 = gp_Ax1(self.rot.Location(), self.rot.YDirection())
        elif axs == "z":
            ax1 = gp_Ax1(self.rot.Location(), self.rot.Direction())
        else:
            ax1 = gp_Ax1(self.rot.Location(), self.rot.Direction())
        rot = self.rot.Rotated(ax1, np.deg2rad(deg))
        trf = gp_Trsf()
        trf.SetDisplacement(self.rot, rot)
        self.axs.Transform(trf)

    def RotateAxis_Ax(self, ax=gp_Ax3(), deg=0.0, axs="z"):
        if axs == "x":
            ax1 = gp_Ax1(ax.Location(), ax.XDirection())
        elif axs == "y":
            ax1 = gp_Ax1(ax.Location(), ax.YDirection())
        elif axs == "z":
            ax1 = gp_Ax1(ax.Location(), ax.Direction())
        else:
            ax1 = gp_Ax1(ax.Location(), ax.Direction())
        rot = ax.Rotated(ax1, np.deg2rad(deg))
        trf = gp_Trsf()
        trf.SetDisplacement(ax, rot)
        ax.Transform(trf)

    def SurfCurvature(self, nxy=[200, 200], lxy=[450, 450], rxy=[700, 0], sxy=[0, 0]):
        px = np.linspace(-1, 1, int(nxy[0])) * lxy[0] / 2
        py = np.linspace(-1, 1, int(nxy[1])) * lxy[1] / 2
        mesh = np.meshgrid(px, py)
        surf_x = curvature(mesh[0], r=rxy[0], s=sxy[0])
        surf_y = curvature(mesh[1], r=rxy[1], s=sxy[1])
        data = surf_x + surf_y
        self.surf, self.surf_pts = surf_spl_pcd(*mesh, data)
        trf = gp_Trsf()
        trf.SetTransformation(self.axs, gp_Ax3())
        self.surf.Location(TopLoc_Location(trf))

    def load_rim(self, rimfile="pln.rim"):
        data = np.loadtxt(rimfile, skiprows=2)
        pts = []
        for xy in data:
            pts.append(gp_Pnt(*xy, 0))
        self.rim = make_polygon(pts, closed=True)

    def load_mat(self, sfcfile="pln_mat.sfc"):
        xs, ys, xe, ye = [float(v) for v in getline(sfcfile, 2).split()]
        nx, ny = [int(v) for v in getline(sfcfile, 3).split()]
        px = np.linspace(xs, xe, nx)
        py = np.linspace(ys, ye, ny)
        mesh = np.meshgrid(px, py)
        data = np.loadtxt(sfcfile, skiprows=3).T
        self.surf, self.surf_pts = surf_spl_pcd(*mesh, data)

    def export_rim_2d(self, rimfile="pln1.rim", name="pln1-rim"):
        rim_2d = dispocc.proj_rim_pln(self, self.rim, self.pln, self.axs)
        fp = open(rimfile, "w")
        fp.write(' {:s}\n'.format(name))
        fp.write('{:12d}{:12d}{:12d}\n'.format(1, 1, 1))
        rim_tmp = gp_Pnt()
        for i, e in enumerate(Topo(rim_2d).edges()):
            e_curve, u0, u1 = BRep_Tool.Curve(e)
            print(i, e, u0, u1)
            if i != 0 and rim_tmp == e_curve.Value(u0):
                u_range = np.linspace(u0, u1, 50)
                rim_tmp = e_curve.Value(u1)
                p = e_curve.Value(u0)
                p.Transform(set_trf(self.axs, gp_Ax3()))
                data = [p.X(), p.Y()]
                print(0, p, u_range[0], u_range[-1])
            elif i != 0 and rim_tmp == e_curve.Value(u1):
                u_range = np.linspace(u1, u0, 50)
                rim_tmp = e_curve.Value(u0)
                p = e_curve.Value(u1)
                p.Transform(set_trf(self.axs, gp_Ax3()))
                data = [p.X(), p.Y()]
                print(1, p, u_range[0], u_range[-1])
            else:
                u_range = np.linspace(u0, u1, 50)
                rim_tmp = e_curve.Value(u1)
                p = e_curve.Value(u0)
                p.Transform(set_trf(self.axs, gp_Ax3()))
                data = [p.X(), p.Y()]
                print(2, p, u_range[0], u_range[-1])
            fp.write(''.join([float_to_string(val) for val in data]) + '\n')
            for u in u_range[1:]:
                p = e_curve.Value(u)
                p.Transform(set_trf(self.axs, gp_Ax3()))
                data = [p.X(), p.Y()]
                fp.write(''.join([float_to_string(val)
                                  for val in data]) + '\n')
        fp.close()

    def export_sfc1_axs(self, sfcfile="pln1_mat.sfc", name="pln1 Mat"):
        surf = BRep_Tool.Surface(self.surf)

        trf = set_trf(self.axs, gp_Ax3())
        xy0 = dispocc.proj_pnt_pln(self, surf.Value(0, 0), self.pln, self.axs)
        xy1 = dispocc.proj_pnt_pln(self, surf.Value(1, 1), self.pln, self.axs)
        xy0.Transform(trf)
        xy1.Transform(trf)

        m2_trf = set_trf(gp_Ax3(), self.axs)
        m2_pln = BRep_Tool.Surface(self.pln)
        for px in np.linspace(-100, 100, 10):
            for py in np.linspace(-100, 100, 10):
                p0 = gp_Pnt(px, py, 0).Transformed(m2_trf)
                p1 = obj.proj_pnt_pln(p0, self.surf, self.axs)

        #ix0, ix1 = m2.surf_pts.LowerRow(), m2.surf_pts.UpperRow()
        #iy0, iy1 = m2.surf_pts.LowerCol(), m2.surf_pts.UpperCol()
        #xy0 = m2.surf_pts.Value(ix0, iy0).Transformed(trf)
        #xy1 = m2.surf_pts.Value(ix1, iy1).Transformed(trf)
        nx, ny = 200, 200
        xs, xe = xy0.X(), xy1.X()
        ys, ye = xy0.Y(), xy1.Y()
        fp = open(sfcfile, "w")
        fp.write(" {} \n".format(name))
        fp.write(" {:.2e} {:.2e} {:.2e} {:.2e}\n".format(xs, ys, xe, ye))
        fp.write(" {:d} {:d}\n".format(nx, ny))
        for ix in np.linspace(0, 1, nx):
            for iy in np.linspace(0, 1, ny):
                p0 = surf.Value(ix, iy)
                p1 = dispocc.proj_pnt_pln(self, p0, self.pln, self.axs)
                pz = p1.Transformed(trf)
                z = p0.Distance(p1)
                fp.write(" {:.5e} ".format(z))
            fp.write("\n")
        fp.close()
        print(xy0)

    def export_sfc2_axs(self, nxy=[200, 200], rx=[-250, 250], ry=[-250, 250], sfcfile="pln1_mat.sfc"):
        trsf = set_trf(gp_Ax3(), self.axs)
        nx, ny = nxy
        xs, xe = rx
        ys, ye = ry
        plnx = np.linspace(xs, xe, nx)
        plny = np.linspace(ys, ye, ny)
        mesh = np.meshgrid(plnx, plny)
        data = np.zeros_like(mesh[0])
        for (ix, iy), x in np.ndenumerate(data):
            px, py = mesh[0][ix, iy], mesh[1][ix, iy]
            p0 = gp_Pnt(px, py, 0).Transformed(trsf)
            p1 = dispocc.proj_pnt_pln(self, p0, self.surf, self.axs)
            z = p0.Distance(p1)
            data[ix, iy] = z
        grasp_sfc(mesh, data, sfcfile)

    def reflect_beam(self, beam0=gp_Ax3(), tr=0):
        v0 = dir_to_vec(beam0.Direction())
        v1 = dir_to_vec(beam0.XDirection())
        surf = BRep_Tool.Surface(self.surf)
        ray = Geom_Line(beam0.Axis())
        uvw = GeomAPI_IntCS(ray, surf).Parameters(1)
        u, v, w = uvw
        p1, vx, vy = gp_Pnt(), gp_Vec(), gp_Vec()
        GeomLProp_SurfaceTool.D1(surf, u, v, p1, vx, vy)
        vz = vx.Crossed(vy)
        if vz.Dot(v0) > 0:
            vz.Reverse()
        vx.Normalize()
        vy.Normalize()
        vz.Normalize()
        self.beam = gp_Ax3(p1,
                           vec_to_dir(v0.Reversed()), vec_to_dir(v1.Reversed()))
        self.norm = gp_Ax3(p1, vec_to_dir(vz), vec_to_dir(vx))
        if tr == 0:
            self.beam.Mirror(self.norm.Ax2())
            if self.beam.Direction().Dot(self.norm.Direction()) < 0:
                self.beam.ZReverse()
        elif tr == 1:
            self.beam.ZReverse()
            # if self.beam.Direction().Dot(self.norm.Direction()) < 0:
            #    self.beam.ZReverse()
        # print(self.beam.Direction().Dot(self.norm.Direction()))


if __name__ == '__main__':
    argvs = sys.argv
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", dest="dir", default="./")
    parser.add_argument("--pxyz", dest="pxyz",
                      default=[0.0, 0.0, 0.0], type="float", nargs=3)
    opt = parser.parse_args()
    print(opt, argvs)
    
    obj = dispocc(touch=True)

    surf0 = SurfOCC()
    obj.display.DisplayShape(surf0.surf, transparency=0.9)

    obj.show_axs_pln()
    obj.show()
