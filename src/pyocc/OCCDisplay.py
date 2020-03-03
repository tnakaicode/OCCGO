import numpy as np
import matplotlib.pyplot as plt
import sys
import time
import os

from .OCCQt import Viewer

from OCC.Display.SimpleGui import init_display
from PyQt5.QtWidgets import QDialog, QLabel, QLineEdit, QPushButton, QVBoxLayout


class OCCDisplay (Viewer):

    def __init__(self):
        self.display, self.start_display, self.add_menu, self.add_function_to_menu = init_display()
        self.display_obj = self.display.View.View()
        super(OCCDisplay, self).__init__()
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
            name = rootname + ".png"
        else:
            name = rootname + ".png"

        print(name)
        self.display.View.Dump(name)
