from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCC.Core.AIS import AIS_Shape
from OCC.Core.Quantity import Quantity_NOC_BLACK
from OCC.Core.gp import gp_XY

from OCC.Display.SimpleGui import init_display

display, start_display, add_menu, add_function_to_menu = init_display()

myBox = BRepPrimAPI_MakeBox(60, 60, 50).Shape()

view = display.View
colorscale = view.ColorScale().GetObject()

aMinRange = colorscale.GetMin()
aMaxRange = colorscale.GetMax()
aNbIntervals = colorscale.GetNumberOfIntervals()
aTextHeight = colorscale.GetTextHeight()
labPosition = colorscale.GetLabelPosition()
position = gp_XY(colorscale.GetXPosition(), colorscale.GetYPosition())
title = colorscale.GetTitle()

view.ColorScaleDisplay()
display.DisplayShape(myBox)
start_display()
