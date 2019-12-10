import numpy as np
import matplotlib.pyplot as plt
import sys, time, os
sys.path.append(os.path.join(".."))

from src.pyocc.OCCQt import Viewer
from src.pyocc.OCCDisplay import OCCDisplay

from OCC.Display.SimpleGui import init_display
from OCC.gp import gp_Pln
from OCC.gp import gp_Pnt, gp_Vec, gp_Ax1, gp_Ax2, gp_Ax3
from OCC.TopoDS  import TopoDS_Compound
from OCC.BOPAlgo import BOPAlgo_MakerVolume, BOPAlgo_Builder
from OCC.BRep    import BRep_Builder
from OCC.BRepAlgoAPI import BRepAlgoAPI_Cut
from OCC.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCC.Graphic3d   import Graphic3d_EF_PDF
from OCCUtils.Topology  import Topo
from OCCUtils.Construct import make_box, make_face
from OCCUtils.Construct import vec_to_dir

class OCCObject (OCCDisplay):

    def __init__(self):
        super(OCCObject, self).__init__()
        print (self.app, self.vi, self.wi)
    
    def Test(self):
        self.display.DisplayShape(make_box(100, 100, 100))
        self.display.DisplayShape(gp_Pnt())
        self.display.FitAll()
        self.start_display()
    
    def Cap (self):
        menu_name = "ScreenCapture"
        self.add_menu (menu_name)
        self.add_function_to_menu (menu_name, self.export_cap)

if __name__ == "__main__":
    obj = OCCObject()
    obj.Cap ()
    obj.Test()