from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCC.Core.AIS import AIS_ColorScale
from OCC.Core.Graphic3d import Graphic3d_ZLayerId_TopOSD, Graphic3d_TMF_2d
from OCC.Core.gp import gp_XY, gp_Pnt

from OCC.Display.SimpleGui import init_display

display, start_display, add_menu, add_function_to_menu = init_display()

myBox = BRepPrimAPI_MakeBox(60, 60, 50).Shape()

colorscale = AIS_ColorScale()

# colorscale properties
aMinRange = colorscale.GetMin()
aMaxRange = colorscale.GetMax()
aNbIntervals = colorscale.GetNumberOfIntervals()
aTextHeight = colorscale.GetTextHeight()
labPosition = colorscale.GetLabelPosition()
position = gp_XY(colorscale.GetXPosition(), colorscale.GetYPosition())
title = colorscale.GetTitle()

# colorscale display
colorscale.SetSize(300, 300)
colorscale.SetRange(0.0, 10.0)
colorscale.SetNumberOfIntervals(10)

colorscale.SetZLayer(Graphic3d_ZLayerId_TopOSD)
colorscale.SetTransformPersistence(Graphic3d_TMF_2d, gp_Pnt(-1, -1, 0))
colorscale.SetToUpdate()

display.Context.Display(colorscale, True)
display.DisplayShape(myBox, update=True)

start_display()
