from src.base import plotocc
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
from linecache import getline, clearcache
from optparse import OptionParser

sys.path.append(os.path.join("./"))
from src.base import plotocc

import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)

from PyQt5.QtWidgets import QApplication, qApp
from PyQt5.QtWidgets import QDialog, QCheckBox
from PyQt5.QtWidgets import QFileDialog
from PyQt5.QtGui import QMovie
from PyQt5.QtCore import QByteArray

from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.TopoDS import TopoDS_Shape, TopoDS_Compound
from OCC.Core.TopLoc import TopLoc_Location
from OCCUtils.Topology import Topo
from OCCUtils.Construct import make_box, make_line, make_wire, make_edge
from OCCUtils.Construct import make_plane, make_polygon
from OCCUtils.Construct import point_to_vector, vector_to_point
from OCCUtils.Construct import dir_to_vec, vec_to_dir

class OCCMovie(plotocc):
    
    def __init__(self, disp=True, touch=True):
        super().__init__(disp=disp, touch=touch)
        self.MovieMenu()
    
    def MovieMenu(self):
        self.add_menu("Movie")
        self.add_function("Movie", self.export_cap)
        self.add_function("Movie", self.movie_start)
        self.add_function("Movie", self.movie_stop)

    
        
    def movie_start(self):
        """sart animnation"""
        # use an animated gif file you have in the working folder
        # or give the full file path
        self.movie = QMovie(self.tempname + ".gif", QByteArray(), self.wi) 
        self.movie.setCacheMode(QMovie.CacheAll) 
        self.movie.setSpeed(100) 
        #self.movie_screen.setMovie(self.movie) 
        self.movie.start()
        
    def movie_stop(self):
        """stop the animation"""
        self.movie.stop()
        
if __name__ == '__main__':
    argvs = sys.argv
    parser = OptionParser()
    parser.add_option("--dir", dest="dir", default="./")
    parser.add_option("--pxyz", dest="pxyz",
                      default=[0.0, 0.0, 0.0], type="float", nargs=3)
    opt, argc = parser.parse_args(argvs)
    print(opt, argc)
    
    # https://www.daniweb.com/programming/software-development/threads/116827/pyqt-and-gif

    obj = OCCMovie(touch=True)
    obj.show_axs_pln()
    obj.show()
