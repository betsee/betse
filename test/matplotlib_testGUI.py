#!/usr/bin/env python3

import sys

from PySide import QtGui
from PySide import QtCore

import matplotlib

matplotlib.use('Qt4Agg')
matplotlib.rcParams['backend.qt4']='PySide'

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib.backend_bases import key_press_handler

import matplotlib.pyplot as plt

import random

class Window(QtGui.QDialog):
    def __init__(self, parent=None):
        super(Window, self).__init__(parent)

        self.figure = plt.figure()
        self.canvas = FigureCanvas(self.figure)


        self.toolbar = NavigationToolbar(self.canvas, self)
        self.toolbar.hide()

        # Just some button
        self.button = QtGui.QPushButton('Plot')
        self.button.clicked.connect(self.plot)

        self.button1 = QtGui.QPushButton('Zoom')
        self.button1.clicked.connect(self.zoom)

        self.button2 = QtGui.QPushButton('Pan')
        self.button2.clicked.connect(self.pan)

        self.button3 = QtGui.QPushButton('Home')
        self.button3.clicked.connect(self.home)


        # set the layout
        layout = QtGui.QVBoxLayout()
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)
        layout.addWidget(self.button)
        layout.addWidget(self.button1)
        layout.addWidget(self.button2)
        layout.addWidget(self.button3)
        self.setLayout(layout)

    def home(self):
        self.toolbar.home()
    def zoom(self):
        self.toolbar.zoom()
    def pan(self):
        self.toolbar.pan()

    def plot(self):
        data = [random.random() for i in range(25)]
        ax = self.figure.add_subplot(111)
        ax.hold(False)
        ax.plot(data, '*-')
        self.canvas.draw()

    # def on_key_press(self, event):
    #     print('you pressed', event.key)
    #     # implement the default mpl key press events described at
    #     # http://matplotlib.org/users/navigation_toolbar.html#navigation-keyboard-shortcuts
    #     key_press_handler(event, self.canvas, self.toolbar)

if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)

    main = Window()
    main.setWindowTitle('Simple QTpy and MatplotLib example with Zoom/Pan')
    main.show()

    sys.exit(app.exec_())