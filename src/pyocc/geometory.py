from OCC.Display.SimpleGui import init_display
from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Ax3
from OCCUtils.Topology import Topo
from OCCUtils.Construct import make_wire, make_edge, make_plane, make_circle, make_line
from OCC.Core.BRep import BRep_Tool
from OCC.Core.TopTools import (TopTools_ListOfShape,
                          TopTools_ListIteratorOfListOfShape,
                          TopTools_IndexedDataMapOfShapeListOfShape)
from OCC.Core.TopoDS import (topods, TopoDS_Wire, TopoDS_Vertex, TopoDS_Edge,
                        TopoDS_Face, TopoDS_Shell, TopoDS_Solid,
                        TopoDS_Compound, TopoDS_CompSolid, topods_Edge,
                        topods_Vertex, TopoDS_Iterator)

from .export import export_STEPFile_single


class OCC_objct (object):
    def __init__(self, display):
        self.display = display
        self.obj = {}

    def display_obj(self, name, size=100):
        self.display_shapes(self.obj[name], "pts")
        self.display.DisplayMessage(self.obj[name]["pnt"], name)
        self.display.DisplayVector(
            self.obj[name]["xyz"][0].Scaled(size), self.obj[name]["pnt"])
        self.display.DisplayVector(
            self.obj[name]["xyz"][1].Scaled(size), self.obj[name]["pnt"])
        self.display.DisplayVector(
            self.obj[name]["xyz"][2].Scaled(size), self.obj[name]["pnt"])
        self.display.DisplayShape(self.obj[name]["frm"])
        self.display.DisplayShape(self.obj[name]["pln"])

        p = self.obj[name]["pnt"]
        v0 = gp_Vec(p.XYZ()) + self.obj[name]["xyz"][0].Scaled(50)
        v1 = gp_Vec(p.XYZ()) + self.obj[name]["xyz"][1].Scaled(100)
        self.display.DisplayShape(make_line(p, gp_Pnt(v0.XYZ())))
        self.display.DisplayShape(make_line(p, gp_Pnt(v1.XYZ())))

    def display_shapes(self, meta, name, color=None):
        for sh in meta[name]:
            self._display_shape(sh, color)

    def _display_shape(self, shape, color=None):
        self.display.DisplayShape(shape, color=color)

    def occ_save(self, names, filename="shape.stp", shp_type="pln"):
        shape = []
        for name in names:
            shape.append(self.obj[name][shp_type])
        export_STEPFile_single(shape, filename)


def set_obj(meta, xyz, dx, dy,
            pln_set=[-200, 200, -200, 200], indx=1):
    vx = gp_Vec(*dx).Normalized()
    vy = gp_Vec(*dy).Normalized()
    vz = vx.Crossed(vy)
    meta["pnt"] = gp_Pnt(*xyz)
    meta["xyz"] = [vx, vy, vz]
    meta["pln"] = make_plane(
        meta["pnt"], meta["xyz"][2],
        *pln_set
    )
    meta["pts"] = shape_to_pts(meta["pln"])
    meta["frm"] = set_wire(meta["pts"])
    meta["idx"] = indx


def shape_to_pts(sh):
    return [BRep_Tool().Pnt(topods_Vertex(v)) for v in Topo(sh).vertices()]


def set_wire(pts):
    edge = []
    for i in range(len(pts)):
        i0, i1 = i, (i+1) % len(pts)
        edge.append(make_edge(pts[i0], pts[i1]))
    return make_wire(edge)
