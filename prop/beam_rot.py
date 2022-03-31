import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
from linecache import getline, clearcache
import argparse

sys.path.append(os.path.join("../"))
from src.base_occ import dispocc

import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)

from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_XYZ
from OCC.Core.gp import gp_Trsf, gp_Quaternion
from OCC.Core.BRep import BRep_Tool
from OCC.Core.GeomAPI import GeomAPI_IntCS, GeomAPI_IntSS
from OCCUtils.Construct import point_to_vector, vector_to_point
from OCCUtils.Construct import dir_to_vec, vec_to_dir

if __name__ == '__main__':
    argvs = sys.argv
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", dest="dir", default="./")
    parser.add_argument("--pxyz", dest="pxyz",
                      default=[0.0, 0.0, 0.0], type=float, nargs=3)
    opt = parser.parse_args()
    print(opt, argvs)

    obj = dispocc(touch=True)
    axs = gp_Ax3(gp_Pnt(100, -100, 200),
                 gp_Dir(0, 0.5, 1.0),
                 gp_Dir(0.5, 1.0, 0))
    qp = gp_Quaternion()
    qp.SetRotation(gp_Vec(0, 0, 100), gp_Vec(50, 0, 0))
    print(qp)
    obj.show_axs_pln(axs, scale=50, name="axs0")

    trf = gp_Trsf()
    trf.SetTransformation(qp, gp_Vec(0, 0, 1))
    ax1 = axs.Transformed(trf)
    obj.show_axs_pln(ax1, scale=50, name="axs1")

    ax2 = obj.prop_axs(ax1, scale=50)
    obj.show_axs_pln(ax2, scale=50, name="axs2")

    ax3 = obj.rot_axis(ax2, xyz="z", deg=45)
    obj.show_axs_pln(ax3, scale=50)
    pln1 = obj.show_plane(ax3, scale=25, trs=0.9, color="RED")
    srf1 = BRep_Tool.Surface(pln1)
    ax3 = obj.rot_axis(ax3, xyz="x", deg=5)
    obj.show_axs_pln(ax3, scale=50)
    pln2 = obj.show_plane(ax3, scale=25, trs=0.9, color="BLUE")
    srf2 = BRep_Tool.Surface(pln2)

    api = GeomAPI_IntSS(srf1, srf2, 0.1E-3)
    obj.display.DisplayShape(api.Line(1))

    obj.show_axs_pln(scale=25)
    obj.show()
