import numpy as np
import matplotlib.pyplot as plt
import sys
import time
import os

from .OCCQt import Viewer

from OCC.Display.SimpleGui import init_display
from OCC.Core.Graphic3d import (Graphic3d_EF_PDF,
                           Graphic3d_EF_SVG,
                           Graphic3d_EF_TEX,
                           Graphic3d_EF_PostScript,
                           Graphic3d_EF_EnhPostScript)

from PyQt5.QtWidgets import QDialog, QLabel, QLineEdit, QPushButton, QVBoxLayout


class OCCDisplay (Viewer):

    def __init__(self):
        self.display, self.start_display, self.add_menu, self.add_function_to_menu = init_display()
        self.display_obj = self.display.View.View().GetObject()
        super(OCCDisplay, self).__init__()
        self.on_select()

    def SubWindow(self):
        self.w = QDialog()
        self.txt = QLineEdit()
        self.layout = QVBoxLayout()
        self.layout.addWidget(self.txt)
        self.w.setLayout(self.layout)
        self.w.exec_()

    def export_cap(self):
        print(os.getcwd())
        sub = self.SubWindow()
        name = self.txt.text()
        indx = name.split(".")[-1]

        self.display.View.Dump(name)
        if indx == "pdf":
            ef_type = Graphic3d_EF_PDF
        elif indx == "ps":
            ef_type = Graphic3d_EF_PostScript
        elif indx == "svg":
            ef_type = Graphic3d_EF_SVG
        else:
            ef_type = Graphic3d_EF_PDF
        print(name)
        #self.display_obj.Export (name, ef_type)
