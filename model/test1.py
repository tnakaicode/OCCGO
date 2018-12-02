import numpy as np
import matplotlib.pyplot as plt
import json, glob
import sys, time, os
from unwrap.unwrap import unwrap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from linecache import getline, clearcache
from scipy.integrate import simps
from scipy.constants import *
from optparse import OptionParser
sys.path.append(os.path.join('..'))

from OCC.Display.SimpleGui import init_display
from OCC.Core.gp import gp_Pln
from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.TopoDS  import TopoDS_Compound
from OCC.Core.BOPAlgo import BOPAlgo_MakerVolume, BOPAlgo_Builder
from OCC.Core.BRep    import BRep_Builder
from OCC.Extend.ShapeFactory import make_box, make_face
from OCCUtils.Topology  import Topo
from OCCUtils.Construct import vec_to_dir

from src.pyocc.OCCQt import Viewer
from src.pyocc.OCCDisplay import OCCDisplay

if __name__ == "__main__":
    print("ok")

    display, start_display, add_menu, add_function_to_menu = init_display()

    display.DisplayShape(make_box(10,10,10))

    display.FitAll()
    start_display()
