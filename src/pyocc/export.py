import numpy as np
import sys, time, os

from OCC.XSControl      import XSControl_Writer, XSControl_WorkSession
from OCC.STEPCAFControl import STEPCAFControl_Writer
from OCC.STEPCAFControl import Handle_STEPCAFControl_ExternFile
from OCC.STEPControl    import STEPControl_Writer
from OCC.STEPControl    import (STEPControl_AsIs,
                                STEPControl_ManifoldSolidBrep,
                                STEPControl_FacetedBrep,
                                STEPControl_ShellBasedSurfaceModel,
                                STEPControl_GeometricCurveSet)
from OCC.Interface   import Interface_Static_SetCVal
from OCC.IFSelect    import IFSelect_RetDone
from OCC.TDF import TDF_LabelSequence, TDF_Label, TDF_Tool
from OCC.TDataStd import Handle_TDataStd_Name, TDataStd_Name_GetID
from OCC.TDataStd import TDataStd_Name
from OCC.TCollection import TCollection_AsciiString
from OCC.TCollection import TCollection_ExtendedString
from OCC.TDocStd import TDocStd_Document, Handle_TDocStd_Document
from OCC.XCAFApp import XCAFApp_Application
from OCC.XCAFDoc import (XCAFDoc_DocumentTool_ShapeTool,
                         XCAFDoc_DocumentTool_ColorTool,
                         XCAFDoc_DocumentTool_LayerTool,
                         XCAFDoc_DocumentTool_MaterialTool)

from OCC.TDF import TDF_Data, TDF_Label, TDF_LabelSequence

class ExportMethod (object):
    
    def __init__(self, tol=1.0E-6):
        self.obj = STEPControl_Writer()
        self.obj.SetTolerance (tol)
        Interface_Static_SetCVal("write.step.schema", "AP214")

        """
        self.obj.PrintStatsTransfer:
        what
        0 gives general statistics (number of translated roots, number of warnings, number of fail messages),
        1 gives root results,
        2 gives statistics for all checked entities,
        3 gives the list of translated entities,
        4 gives warning and fail messages,
        5 gives fail messages only. The use of mode depends on the value of what. If what is 0, mode is ignored. 
        
        If what is 1, 2 or 3, mode defines the following:
        0 lists the numbers of IGES or STEP entities in the respective model
        1 gives the number, identifier, type and result type for each IGES or STEP entity and/or its status (fail, warning, etc.)
        2 gives maximum information for each IGES or STEP entity (i.e. checks)
        3 gives the number of entities per type of IGES or STEP entity
        4 gives the number of IGES or STEP entities per result type and/or status
        5 gives the number of pairs (IGES or STEP or result type and status)
        6 gives the number of pairs (IGES or STEP or result type and status) AND the list of entity numbers in the IGES or STEP model. 
        
        If what is 4 or 5, mode defines the warning and fail messages as follows:
        if mode is 0 all warnings and checks per entity are returned
        if mode is 2 the list of entities per warning is returned. If mode is not set, only the list of all entities per warning is given.
        """
    
    def add_shpe (self, shape):
        """
        STEPControl_AsIs                   translates an Open CASCADE shape to its highest possible STEP representation.
        STEPControl_ManifoldSolidBrep      translates an Open CASCADE shape to a STEP manifold_solid_brep or brep_with_voids entity.
        STEPControl_FacetedBrep            translates an Open CASCADE shape into a STEP faceted_brep entity.
        STEPControl_ShellBasedSurfaceModel translates an Open CASCADE shape into a STEP shell_based_surface_model entity.
        STEPControl_GeometricCurveSet      translates an Open CASCADE shape into a STEP geometric_curve_set entity.
        """
        self.obj.Transfer(shape, STEPControl_AsIs)
    
    def fileout (self, filename):
        status = self.obj.Write(filename)
        assert(status == IFSelect_RetDone)


class ExportCAFMethod (object):
    
    def __init__(self, name="name", tol=1.0E-10):
        self.name = name
        self.step = STEPCAFControl_Writer()
        self.step.SetNameMode(True)
        self.h_doc = Handle_TDocStd_Document()
        self.x_app = XCAFApp_Application.GetApplication().GetObject()
        self.x_app.NewDocument(TCollection_ExtendedString("MDTV-CAF"), self.h_doc)
        self.doc   = self.h_doc.GetObject()
        self.h_shape_tool = XCAFDoc_DocumentTool_ShapeTool(self.doc.Main())
        self.shape_tool   = self.h_shape_tool.GetObject()
        Interface_Static_SetCVal("write.step.schema", "AP214")
        
    def Add (self, shape, name="name"):
        """
        STEPControl_AsIs                   translates an Open CASCADE shape to its highest possible STEP representation.
        STEPControl_ManifoldSolidBrep      translates an Open CASCADE shape to a STEP manifold_solid_brep or brep_with_voids entity.
        STEPControl_FacetedBrep            translates an Open CASCADE shape into a STEP faceted_brep entity.
        STEPControl_ShellBasedSurfaceModel translates an Open CASCADE shape into a STEP shell_based_surface_model entity.
        STEPControl_GeometricCurveSet      translates an Open CASCADE shape into a STEP geometric_curve_set entity.
        """
        label = self.shape_tool.AddShape(shape)
        self.step.Transfer(self.h_doc, STEPControl_AsIs)
    
    def Write (self, filename=None):
        if not filename:
            filename = self.name
        path, ext = os.path.splitext(filename)
        if not ext:
            ext = ".stp"
        status = self.step.Write(path + ext)
        assert(status == IFSelect_RetDone)
        
def export_STEPFile_single(shape, filename, tol=1.0E-6):
    """
    Exports a .stp file containing the input shapes

    Parameters
    ----------
    shape : TopoDS_Shape
    
    filename : string
        The output filename
    """

    step = STEPCAFControl_Writer()
    step.SetNameMode(True)
    step.SetPropsMode(True)
    h_doc = Handle_TDocStd_Document()
    x_app = XCAFApp_Application.GetApplication().GetObject()
    x_app.NewDocument(TCollection_ExtendedString("MDTV-CAF"), h_doc)
    doc   = h_doc.GetObject()
    h_shape_tool = XCAFDoc_DocumentTool_ShapeTool(doc.Main())
    shape_tool   = h_shape_tool.GetObject()
    Interface_Static_SetCVal("write.step.schema", "AP214")
    
    # transfer shapes
    print (filename)
    shape_tool.AddShape(shape)
    step.Transfer(h_doc, STEPControl_AsIs)
    status = step.Write(filename)
    assert(status == IFSelect_RetDone)

def export_STEPFile_name (shapes, filename):
    print (filename)
    step_writer = STEPControl_Writer()
    ws = step_writer.WS().GetObject()
    print (ws.ModeWriteShape())

def write_step(shape, file_name):
    step_writer = STEPControl_Writer()
    step_writer.Transfer(shape, STEPControl_AsIs)
    step_writer.Write(file_name)
