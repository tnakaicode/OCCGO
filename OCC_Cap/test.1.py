import sys
from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from PyQt5.QtCore import *


class MainWindow(QWidget):
    def __init__(self, parent=None):
        super(MainWindow, self).__init__(parent)

        makeWindowButton = QPushButton("&make window")
        makeWindowButton.clicked.connect(self.makeWindow)

        self.label = QLabel()
        self.label.setText("default:")

        layout = QHBoxLayout()
        layout.addWidget(makeWindowButton)
        layout.addWidget(self.label)
        self.setLayout(layout)

    def makeWindow(self):
        subWindow = SubWindow(self)
        subWindow.show()

    # サブウィンドウから実行
    def setParam(self, param):
        self.label.setText(param)

class SubWindow:
    def __init__(self, parent=None):
        self.w = QDialog(parent)
        self.parent = parent

        label = QLabel()
        label.setText('Sub Window')

        self.edit = QLineEdit()

        button = QPushButton('送信')
        button.clicked.connect(self.setParamOriginal)

        layout = QVBoxLayout()
        layout.addWidget(label)
        layout.addWidget(self.edit)
        layout.addWidget(button)

        self.w.setLayout(layout)

    # ここで親ウィンドウに値を渡している
    def setParamOriginal(self):
        self.parent.setParam(self.edit.text())

    def show(self):
        self.w.exec_()

if __name__ == '__main__':
    import sys
    app = QApplication(sys.argv)
    main_window = MainWindow()

    main_window.show()
    sys.exit(app.exec_())
