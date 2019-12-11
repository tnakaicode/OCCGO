import numpy as np
import matplotlib.pyplot as plt
import json
import glob
import sys
import time
import os
import scipy.constants as cnt
import random
from unwrap.unwrap import unwrap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from linecache import getline, clearcache
from scipy.integrate import simps
from optparse import OptionParser

from OCC.Display.SimpleGui import init_display
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Pln, gp_Trsf, gp_Lin
from OCC.Core.gp import gp_Elips, gp_Elips2d
from OCC.Core.gp import gp_Mat
from OCC.Core.BRep import BRep_Tool_PolygonOnTriangulation
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Algo
from OCC.Core.Graphic3d import Graphic3d_ArrayOfPoints
from OCC.Core.AIS import AIS_PointCloud
from OCC.Core.Bnd import Bnd_Box, Bnd_Array1OfSphere


if __name__ == "__main__":
    argvs = sys.argv
    parser = OptionParser()
    opt, argc = parser.parse_args(argvs)
    print(argc, opt)

    display, start_display, add_menu, add_function_to_menu = init_display()

    display.DisplayShape(gp_Pnt())

    pnt_3d = Graphic3d_ArrayOfPoints(100)
    for idx in range(100):
        x = random.uniform(-0.5, 0.5)
        y = random.uniform(-0.5, 0.5)
        z = random.uniform(-0.5, 0.5)
        r = np.sqrt(x**2 + y**2 + z**2)
        ratio = random.uniform(0, 100)
        pnt_3d.AddVertex(x/r*ratio, y/r*ratio, z/r*ratio)
    
    point_cloud = AIS_PointCloud()
    point_cloud.SetPoints(pnt_3d)

    #display.DisplayShape(point_cloud.GetBoundingBox())

    display.FitAll()
    start_display()
