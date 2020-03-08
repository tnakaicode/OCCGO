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

import os
import sys
import time

from OCC.Display.SimpleGui import init_display
from OCC.Core.BRep import BRep_Tool
from OCC.Core.BRepAdaptor import BRepAdaptor_HCurve
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakePolygon
from OCC.Core.BRepFill import BRepFill_CurveConstraint
from OCC.Core.GeomAbs import GeomAbs_C0
from OCC.Core.GeomLProp import GeomLProp_SLProps
from OCC.Core.GeomPlate import (GeomPlate_BuildPlateSurface, GeomPlate_PointConstraint,
                                GeomPlate_MakeApprox)
from OCC.Core.ShapeAnalysis import ShapeAnalysis_Surface
from OCC.Core.gp import gp_Pnt
from OCC.Core.BRepFill import BRepFill_Filling
from OCC.Extend.TopologyUtils import TopologyExplorer, WireExplorer
from OCC.Extend.DataExchange import read_iges_file
from OCC.Extend.ShapeFactory import make_edge, make_n_sided
from OCCUtils.Construct import make_plane


def geom_plate():
    p1 = gp_Pnt(0, 0, 0)
    p2 = gp_Pnt(10, 0, 0)
    p3 = gp_Pnt(10, 0, 9)
    p4 = gp_Pnt(0, 1, 1)
    p5 = gp_Pnt(5, 5, 0)
    edge = [p1, p2, p3, p4]
    edge.append(edge[0])
    poly = [make_edge(edge[i], edge[i + 1]) for i in range(len(edge) - 1)]
    face = make_n_sided(poly)
    return face


if __name__ == "__main__":
    display, start_display, add_menu, add_function_to_menu = init_display()
    display.DisplayShape(geom_plate())
    display.DisplayShape(make_plane(
        extent_x_max=10, extent_x_min=-10, extent_y_min=-10, extent_y_max=10), transparency=0.9)

    display.FitAll()
    start_display()
