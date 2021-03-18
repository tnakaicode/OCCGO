import sys

# too lazy to keep track of QtCore or QtGui
from PyQt5.QtGui import QMovie
from PyQt5.QtCore import QByteArray
from PyQt5.QtWidgets import QApplication, QLabel, QPushButton, QSizePolicy, QVBoxLayout, qApp, QMainWindow
from PyQt5.QtWidgets import QDialog, QCheckBox
from PyQt5.QtWidgets import QFileDialog


class MoviePlayer(QMainWindow):

    def __init__(self, *args):
        QMainWindow.__init__(self, *args)
        # setGeometry(x_pos, y_pos, width, height)
        self.setGeometry(200, 200, 400, 400)
        self.setWindowTitle("QMovie to show animated gif")

        # set up the movie screen on a label

        self.movie_screen = QLabel()
        # expand and center the label
        self.movie_screen.setSizePolicy(QSizePolicy.Expanding,
                                        QSizePolicy.Expanding)

        #btn_start = QPushButton("Start Animation", self)
        # btn_start.clicked.connect(self.start)

        btn_stop = QPushButton("Stop Animation", self)
        btn_stop.clicked.connect(self.stop)

        main_layout = QVBoxLayout()
        main_layout.addWidget(self.movie_screen)
        # main_layout.addWidget(btn_start)
        main_layout.addWidget(btn_stop)
        self.setLayout(main_layout)

        # use an animated gif file you have in the working folder
        # or give the full file path
        self.movie = QMovie("./AG_Dog.gif", QByteArray(), self)
        self.movie.setCacheMode(QMovie.CacheAll)
        self.movie.setSpeed(100)
        self.movie_screen.setMovie(self.movie)
        self.movie.start()

    def start(self):
        """sart animnation"""
        self.movie.start()

    def stop(self):
        """stop the animation"""
        self.movie.stop()
        print(self.movie.fileName())
        self.movie.start()


app = QApplication(sys.argv)
player = MoviePlayer()
player.show()
sys.exit(app.exec_())
