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

sys.path.append(os.path.join("../"))
from src.base import plotocc

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

    obj = plotocc()
    obj.display.DisplayShape(gp_Pnt())

    num = 100000
    pnt_3d = Graphic3d_ArrayOfPoints(num)
    for idx in range(num):
        x = random.uniform(-0.5, 0.5)
        y = random.uniform(-0.5, 0.5)
        z = random.uniform(-0.5, 0.5)
        r = np.sqrt(x**2 + y**2 + z**2)
        ratio = random.uniform(0, 100)
        pnt_3d.AddVertex(x / r * ratio, y / r * ratio, z / r * ratio)

    point_cloud = AIS_PointCloud()
    point_cloud.SetPoints(pnt_3d)
    box = point_cloud.GetBoundingBox()
    print(box)

    ais_context = obj.display.GetContext()
    ais_context.Display(point_cloud, True)
    obj.show()
