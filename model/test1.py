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
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/..")

from OCC.Display.SimpleGui import init_display
from OCC.gp import gp_Pln
from OCC.gp import gp_Pnt, gp_Vec, gp_Ax1, gp_Ax2, gp_Ax3
from OCC.TopoDS  import TopoDS_Compound
from OCC.BOPAlgo import BOPAlgo_MakerVolume, BOPAlgo_Builder
from OCC.BRep    import BRep_Builder
from OCCUtils.Topology  import Topo
from OCCUtils.Construct import make_box, make_face
from OCCUtils.Construct import vec_to_dir

from src.pyocc.OCCQt import Viewer
from src.pyocc.OCCDisplay import OCCDisplay

if __name__ == "__main__":
    print("ok")
