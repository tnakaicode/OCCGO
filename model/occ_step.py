import numpy as np
import sys
import time
import os

from OCC.Display.SimpleGui import init_display
from OCC.Core.gp import gp_Pnt
from OCC.Core.XSControl import XSControl_Writer, XSControl_WorkSession
from OCC.Core.XCAFApp import XCAFApp_Application
from OCC.Core.XCAFDoc import (XCAFDoc_DocumentTool_ShapeTool,
                              XCAFDoc_DocumentTool_ColorTool,
                              XCAFDoc_DocumentTool_LayerTool,
                              XCAFDoc_DocumentTool_MaterialTool)
from OCC.Core.STEPCAFControl import STEPCAFControl_Writer
from OCC.Core.STEPConstruct import stepconstruct_FindEntity
from OCC.Core.STEPControl import STEPControl_Writer
from OCC.Core.STEPControl import (STEPControl_AsIs,
                                  STEPControl_ManifoldSolidBrep,
                                  STEPControl_FacetedBrep,
                                  STEPControl_ShellBasedSurfaceModel,
                                  STEPControl_GeometricCurveSet)
from OCC.Core.Interface import Interface_Static_SetCVal
from OCC.Core.IFSelect import IFSelect_RetDone
from OCC.Core.TDF import TDF_LabelSequence, TDF_Label, TDF_Tool, TDF_Data
from OCC.Core.TDataStd import TDataStd_Name, TDataStd_Name_GetID
from OCC.Core.TCollection import TCollection_AsciiString, TCollection_HAsciiString
from OCC.Core.TCollection import TCollection_ExtendedString
from OCC.Core.TDocStd import TDocStd_Document
from OCCUtils.Construct import make_box
from OCCUtils.Construct import make_plane, make_vertex, make_circle


class ExportCAFMethod (object):

    def __init__(self, name="name", tol=1.0E-10):
        self.name = name
        self.schema = 'AP214'
        self.assembly_mode = 1

        self.stp = STEPControl_Writer()
        self.stp.SetTolerance(tol)
        self.app = self.stp.WS().TransferWriter().FinderProcess()
        Interface_Static_SetCVal('write.step.schema', self.schema)
        Interface_Static_SetCVal('write.step.unit', 'MM')
        Interface_Static_SetCVal(
            'write.step.assembly', str(self.assembly_mode))

        # Interface_Static_SetCVal ("write.step.schema","AP203")
        # Interface_Static_SetIVal ("write.step.schema", 3)
        # Interface_Static_SetRVal

    def Add(self, shape, name="name"):
        """
        STEPControl_AsIs                   translates an Open CASCADE shape to its highest possible STEP representation.
        STEPControl_ManifoldSolidBrep      translates an Open CASCADE shape to a STEP manifold_solid_brep or brep_with_voids entity.
        STEPControl_FacetedBrep            translates an Open CASCADE shape into a STEP faceted_brep entity.
        STEPControl_ShellBasedSurfaceModel translates an Open CASCADE shape into a STEP shell_based_surface_model entity.
        STEPControl_GeometricCurveSet      translates an Open CASCADE shape into a STEP geometric_curve_set entity.
        """
        Interface_Static_SetCVal('write.step.product.name', name)
        self.stp.Transfer(shape, STEPControl_AsIs)
        item = stepconstruct_FindEntity(self.app, shape)
        item.SetName(TCollection_HAsciiString(name))

    def Write(self, filename=None):
        if not filename:
            filename = self.name
        path, ext = os.path.splitext(filename)
        if not ext:
            ext = ".stp"
        status = self.stp.Write(path + ext)
        assert(status == IFSelect_RetDone)


if __name__ == "__main__":
    display, start_display, add_menu, add_function_to_menu = init_display()
    display.DisplayShape(gp_Pnt())
    root = ExportCAFMethod(name="root")
    root.Add(make_vertex(gp_Pnt()), name="pnt")
    root.Add(make_plane(center=gp_Pnt(0, 0, 0)), name="pln0")
    root.Add(make_plane(center=gp_Pnt(0, 0, 100)), name="pln1")
    root.Add(make_plane(center=gp_Pnt(0, 0, 200)), name="pln2")
    root.Add(make_circle(gp_Pnt(0, 0, 0), 100), name="circle")
    root.Add(make_box(100, 100, 100), name="box001")
    root.Write()

    display.FitAll()
    # start_display()
