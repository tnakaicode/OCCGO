import numpy as np
import matplotlib.pyplot as plt

from OCC.Display.SimpleGui import init_display
from OCC.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.gp import gp_Pln, gp_Trsf, gp_Lin

from OCC.Extend.ShapeFactory import make_box
# from OCCUtils.Construct import make_box

if __name__ == "__main__":
    display, start_display, add_menu, add_function_to_menu = init_display()

    pnt = gp_Pnt()
    display.DisplayShape(pnt)
    display.DisplayShape(make_box(pnt, 100, 100, 100), transparency=0.1)
    display.DisplayShape(make_box(gp_Pnt(50,50,50), 10, 10, 10))

    display.FitAll()
    start_display()
