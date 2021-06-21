import numpy as np
import matplotlib.pyplot as plt
import sys
import time
import os

from .OCCQt import Viewer
from ..base_occ import dispocc

from OCC.Display.SimpleGui import init_display
from PyQt5.QtWidgets import QDialog, QLabel, QLineEdit, QPushButton, QVBoxLayout


class OCCDisplay (Viewer, dispocc):

    def __init__(self):
        dispocc.__init__(self)
        self.display_obj = self.display.View.View()
        Viewer.__init__(self)
        self.on_select()

    def SubWindow(self):
        self.w = QDialog()
        self.txt = QLineEdit()
        button = QPushButton("enter")
        button.clicked.connect(self.w.close)
        button.resize(button.sizeHint())
        button.resize(200, 32)
        button.move(80, 60)
        layout = QVBoxLayout()
        layout.addWidget(self.txt)
        layout.addWidget(button)
        self.w.setLayout(layout)
        self.w.exec_()

    def export_cap(self):
        print(os.getcwd())
        self.SubWindow()
        rootname, ext_name = os.path.splitext(self.txt.text())
        if ext_name == "":
            name = self.tmpdir + rootname + ".png"
        else:
            name = self.tmpdir + rootname + ".png"

        print(name)
        self.display.View.Dump(name)
