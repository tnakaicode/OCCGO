# Copyright 2009-2016 Jelle Feringa (jelleferinga@gmail.com)
##
# This file is part of pythonOCC.
##
# pythonOCC is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
##
# pythonOCC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
##
# You should have received a copy of the GNU Lesser General Public License
# along with pythonOCC.  If not, see <http://www.gnu.org/licenses/>.
from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import time

sys.path.append(os.path.join("../"))
from src.base import plotocc, PlotBase

from OCC.Display.SimpleGui import init_display
from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.BRep import BRep_Tool
from OCC.Core.BRepAdaptor import BRepAdaptor_HCurve
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakePolygon
from OCC.Core.BRepFill import BRepFill_Filling, BRepFill_CurveConstraint
from OCC.Core.GeomAbs import GeomAbs_C0
from OCC.Core.GeomLProp import GeomLProp_SLProps
from OCC.Core.ShapeAnalysis import ShapeAnalysis_Surface
from OCC.Core.TColgp import TColgp_Array2OfPnt, TColgp_Array1OfPnt
from OCC.Core.GeomAPI import GeomAPI_PointsToBSplineSurface
from OCC.Core.GeomAPI import GeomAPI_PointsToBSpline
from OCC.Core.GeomAPI import GeomAPI_ProjectPointOnSurf
from OCC.Core.GeomAbs import GeomAbs_C2, GeomAbs_C0, GeomAbs_G1, GeomAbs_G2
from OCC.Core.BRep import BRep_Tool_Surface, BRep_Builder
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeFace
from OCC.Core.BRepExtrema import BRepExtrema_DistShapeShape

from OCC.Extend.TopologyUtils import TopologyExplorer, WireExplorer
from OCC.Extend.ShapeFactory import make_face, make_vertex
from OCC.Extend.DataExchange import read_iges_file
from OCC.Extend.ShapeFactory import make_edge, make_n_sided


def curvature(px, r, s):
    """( x + sx )**2 / 2*rx + ( y + sy )**2 / 2*ry"""
    if (r == 0):
        py = np.zeros_like(px + s)
    else:
        py = (px + s)**2 / (2 * r)
    return py


def geom_plate():
    p1 = gp_Pnt(0, 0, 0)
    p2 = gp_Pnt(0, 10, 0)
    p3 = gp_Pnt(0, 10, 9)
    p4 = gp_Pnt(1, 0, 11)
    p5 = gp_Pnt(0, 5, 5)
    edge = [p1, p2, p3, p4]
    edge.append(edge[0])
    poly = [make_edge(edge[i], edge[i + 1]) for i in range(len(edge) - 1)]
    face = make_n_sided(poly)
    return face


class OCC_Display (object):

    def __init__(self):
        self.display, self.start_display, self.add_menu, self.add_function_to_menu = init_display()


class Surf (plotocc):

    def __init__(self, disp=True):
        plotocc.__init__(self, disp=disp, touch=True)
        self.nx = 10
        self.ny = 20
        self.lx = 150
        self.ly = 200
        self.ini_mesh()

    def ini_mesh(self):
        px = np.linspace(-1, 1, self.nx) * self.lx / 2
        py = np.linspace(-1, 1, self.ny) * self.ly / 2
        self.mesh = np.meshgrid(px, py)

    def gen_surf(self, rxy=[100, 200], sxy=[0, 0]):
        surf_x = curvature(self.mesh[0], rxy[0], sxy[0])
        surf_y = curvature(self.mesh[1], rxy[1], sxy[1])
        self.surf = surf_x + surf_y

        self.pnt_2d = TColgp_Array2OfPnt(1, self.nx, 1, self.ny)
        for idx_row in range(self.pnt_2d.LowerRow(), self.pnt_2d.UpperRow() + 1):
            for idx_col in range(self.pnt_2d.LowerCol(), self.pnt_2d.UpperCol() + 1):
                row = idx_row - 1
                col = idx_col - 1
                pnt = gp_Pnt(
                    self.mesh[0][col, row],
                    self.mesh[1][col, row],
                    self.surf[col, row]
                )
                self.pnt_2d.SetValue(idx_row, idx_col, pnt)
                # self.display.DisplayShape(pnt)

    def gen_face(self):
        curv = GeomAPI_PointsToBSplineSurface(
            self.pnt_2d, 3, 8, GeomAbs_G2, 0.001).Surface()
        surf = BRepBuilderAPI_MakeFace(curv, 1e-6).Face()
        self.display.DisplayShape(surf)

    def gen_fill(self):
        # curv = GeomAPI_PointsToBSplineSurface(
        #    self.pnt_2d, 3, 8, GeomAbs_G2, 0.001).Surface()
        # surf = BRepBuilderAPI_MakeFace(curv, 1e-6).Face()
        # self.display.DisplayShape(surf)

        for idx_row in range(self.pnt_2d.LowerRow(), self.pnt_2d.UpperRow()):
            for idx_col in range(self.pnt_2d.LowerCol(), self.pnt_2d.UpperCol()):
                i00, i01 = idx_row, idx_row + 1
                i10, i11 = idx_col, idx_col + 1
                p1 = self.pnt_2d.Value(i00, i10)
                p2 = self.pnt_2d.Value(i01, i10)
                p3 = self.pnt_2d.Value(i01, i11)
                p4 = self.pnt_2d.Value(i00, i11)
                edge = [p1, p2, p3, p4]
                edge.append(edge[0])
                poly = [make_edge(edge[i], edge[i + 1])
                        for i in range(len(edge) - 1)]
                face = make_n_sided(poly)
                self.display.DisplayShape(face)


if __name__ == "__main__":
    print("ok")
    obj = Surf(disp=True)
    obj.gen_surf(sxy=[20, 10])
    obj.gen_fill()
    obj.show()
