##Copyright 2010-2014 Thomas Paviot (tpaviot@gmail.com)
##
##This file is part of pythonOCC.
##
##pythonOCC is free software: you can redistribute it and/or modify
##it under the terms of the GNU Lesser General Public License as published by
##the Free Software Foundation, either version 3 of the License, or
##(at your option) any later version.
##
##pythonOCC is distributed in the hope that it will be useful,
##but WITHOUT ANY WARRANTY; without even the implied warranty of
##MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##GNU Lesser General Public License for more details.
##
##You should have received a copy of the GNU Lesser General Public License
##along with pythonOCC.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import print_function

import sys

from OCC.STEPControl import STEPControl_Reader
from OCC.IGESControl import IGESControl_Reader
from OCC.IFSelect import IFSelect_RetDone, IFSelect_ItemsByEntity
from OCC.Quantity import Quantity_Color, Quantity_TOC_RGB
from OCC.TDocStd import Handle_TDocStd_Document
from OCC.XCAFApp import XCAFApp_Application
from OCC.XCAFDoc import (XCAFDoc_DocumentTool_ShapeTool,
                         XCAFDoc_DocumentTool_ColorTool,
                         XCAFDoc_DocumentTool_LayerTool,
                         XCAFDoc_DocumentTool_MaterialTool)
from OCC.STEPCAFControl import STEPCAFControl_Reader
from OCC.TDF import TDF_LabelSequence
from OCC.TCollection import TCollection_ExtendedString

def read_step_file(filename):
    """ read the STEP file and returns a compound
    """
    step_reader = STEPControl_Reader()
    status = step_reader.ReadFile(filename)

    if status == IFSelect_RetDone:  # check status
        failsonly = False
        step_reader.PrintCheckLoad(failsonly, IFSelect_ItemsByEntity)
        step_reader.PrintCheckTransfer(failsonly, IFSelect_ItemsByEntity)

        ok = step_reader.TransferRoot()
        _nbs = step_reader.NbShapes()
        aResShape = step_reader.OneShape()
        print (_nbs)
    else:
        print("Error: can't read file.")
        sys.exit(0)
    return aResShape

def read_step_file_shapes(filename):
    _shapes = []
    
    # create an handle to a document
    h_doc = Handle_TDocStd_Document()

    # Create the application
    app = XCAFApp_Application.GetApplication().GetObject()
    app.NewDocument(TCollection_ExtendedString("MDTV-CAF"), h_doc)
    
    # Get root assembly
    doc = h_doc.GetObject()
    h_shape_tool = XCAFDoc_DocumentTool_ShapeTool(doc.Main())
    
    step_reader = STEPCAFControl_Reader()
    step_reader.SetNameMode(True)
    
    status = step_reader.ReadFile(filename)
    if status == IFSelect_RetDone:
        step_reader.Transfer(doc.GetHandle())
    
    labels = TDF_LabelSequence()
    shape_tool = h_shape_tool.GetObject()
    h_shape_tool.GetObject().GetFreeShapes(labels)
    
    print("Number of shapes at root :%i" % labels.Length())
    for i in range(labels.Length()):
        label = labels.Value(i+1)
        a_shape = h_shape_tool.GetObject().GetShape(label)
        _shapes.append(a_shape)
    return _shapes

def read_step (file_name):
    step_reader = STEPControl_Reader()
    status = step_reader.ReadFile(file_name)

    if status == IFSelect_RetDone:  # check status
        failsonly = False
        step_reader.PrintCheckLoad(failsonly, IFSelect_ItemsByEntity)
        step_reader.PrintCheckTransfer(failsonly, IFSelect_ItemsByEntity)

        ok = step_reader.TransferRoot(1)
        _nbs = step_reader.NbShapes()
        aResShape = step_reader.Shape(1)
        return aResShape
    else:
        print("Error: can't read file.")
        sys.exit(0)

def read_iges (file_name):
    iges_reader = IGESControl_Reader()
    status = iges_reader.ReadFile(file_name)
    
    if status == IFSelect_RetDone:
        failsonly = False
        iges_reader.PrintCheckLoad(failsonly, IFSelect_ItemsByEntity)
        iges_reader.PrintCheckTransfer(failsonly, IFSelect_ItemsByEntity)
        ok = iges_reader.TransferRoots()
        IgesShape = iges_reader.OneShape()
        return IgesShape
    else:
        print("Error: can't read file.")
        sys.exit(0)
