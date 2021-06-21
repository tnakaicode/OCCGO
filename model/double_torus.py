from OCC.Core.Geom import Geom_ToroidalSurface
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
from linecache import getline, clearcache
from optparse import OptionParser

sys.path.append(os.path.join("./"))
from src.base import plotocc, plot2d

import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)

from OCC.Core.gp import gp_Pnt, gp_Torus, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.TopoDS import TopoDS_Shape, TopoDS_Compound
from OCC.Core.TopLoc import TopLoc_Location
from OCCUtils.Topology import Topo
from OCCUtils.Construct import make_box, make_face, make_line, make_wire, make_edge
from OCCUtils.Construct import make_plane, make_polygon
from OCCUtils.Construct import point_to_vector, vector_to_point
from OCCUtils.Construct import dir_to_vec, vec_to_dir

if __name__ == '__main__':
    argvs = sys.argv
    parser = OptionParser()
    parser.add_option("--dir", dest="dir", default="./")
    parser.add_option("--pxyz", dest="pxyz",
                      default=[0.0, 0.0, 0.0], type="float", nargs=3)
    opt, argc = parser.parse_args(argvs)
    print(opt, argc)

    px = np.linspace(-1, 1, 100) * 100 + 50
    py = np.linspace(-1, 1, 200) * 100 - 50
    mesh = np.meshgrid(px, py)

    p2d = plot2d(aspect="auto")
    p2d.contourf_sub2(mesh, mesh[0], pngname=p2d.tempname + "_plot2d")

    obj = plotocc(touch=True)

    surf1 = make_face(gp_Torus(gp_Ax3(), 1000, 500),
                      0, 2 * np.pi, -np.pi, np.pi)

    surf2 = make_face(gp_Torus(gp_Ax3(), 1500, 750),
                      0, 2 * np.pi, np.pi * 1 / 20, np.pi * 21 / 20)

    surf3 = make_face(gp_Torus(gp_Ax3(), 2000, 1000),
                      0, 2 * np.pi, -np.pi * 15 / 20, -np.pi * 7 / 20)

    obj.display.DisplayShape(surf1, transparency=0.9)
    obj.display.DisplayShape(surf2, transparency=0.9)
    obj.display.DisplayShape(surf3, transparency=0.9)

    obj.show_axs_pln()
    obj.show()
