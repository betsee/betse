#!/usr/bin/env python3

# Import the main necessities for creating GUIs with namespaces
import sys
from PySide import QtCore
from PySide import QtGui

# This creates the main class that is the GUI application

class Example(QtGui.QWidget):

    def __init__(self):      # The self and parent constructor calls....
        super().__init__()

        self.initUi()    # a call to the GUI initialization method.


    def initUi(self):    # this is where the GUI settings are created

        lcd = QtGui.QLCDNumber(self)
        sld = QtGui.QSlider(QtCore.Qt.Vertical, self)

        hbox = QtGui.QHBoxLayout()
        hbox.addWidget(lcd)
        hbox.addWidget(sld)

        self.setLayout(hbox)
        sld.valueChanged.connect(lcd.display)

        # icon = QtGui.QPixmap("plant.png")   # this is how to load an image for a label
        #
        #
        # names = ['clr','back',' ', 'Close', '7', '8', '9', '/', '4','5','6','*','1','2','3','-','0','.','=','+']
        # grid = QtGui.QGridLayout()
        #
        # j = 0
        # pos = [(0, 0), (0, 1), (0, 2), (0, 3),
        #         (1, 0), (1, 1), (1, 2), (1, 3),
        #         (2, 0), (2, 1), (2, 2), (2, 3),
        #         (3, 0), (3, 1), (3, 2), (3, 3 ),
        #         (4, 0), (4, 1), (4, 2), (4, 3)]
        #
        # for i in names:
        #     button = QtGui.QPushButton(i)
        #     if j == 2:
        #         hello = QtGui.QLabel(self)
        #         hello.setPixmap(icon)
        #         grid.addWidget(hello, 0, 2)
        #     else:
        #         grid.addWidget(button, pos[j][0], pos[j][1])
        #     j = j + 1
        #
        # self.setLayout(grid)

        # okBtn = QtGui.QPushButton("OK", self)
        # okBtn.setIcon(QtGui.QIcon("owl.svg"))
        #
        # canBtn = QtGui.QPushButton("Cancel", self)
        # canBtn.setIcon(QtGui.QIcon("ant.svg"))
        #
        # hbox = QtGui.QHBoxLayout()
        # hbox.addStretch(1)
        # hbox.addWidget(okBtn)
        # hbox.addWidget(canBtn)
        #
        # vbox = QtGui.QVBoxLayout()
        # vbox.addStretch(1)
        # vbox.addLayout(hbox)
        #
        # self.setLayout(vbox)

        # self.statusBar()                    # this creates a status bar that continually updates

        self.setGeometry(300,300,500,500)   # this sets the position and size of the application
        self.setWindowTitle("Yolo")         # this titles the main window

        # exitAction = QtGui.QAction(QtGui.QIcon('darth.svg'),'Exit', self)    # this defines an action with an icon
        # exitAction.setShortcut('Ctrl+Q')                # this puts on a keyboard shortcut to the action
        # exitAction.setStatusTip('Exit application')     # this puts a status tip to the action
        # exitAction.triggered.connect(self.close)        # this connects the execution of the action to closing
        #
        # menubar = self.menuBar()                         # this creates a menu bar
        # fileMenu = menubar.addMenu('&File')              # this adds an item to the menubar
        # fileMenu.addAction(exitAction)                   # this adds an action under the item just added
        # fileMenu = menubar.addMenu('&Settings')          # this appends another item in the menubar
        #
        # self.toolbar = self.addToolBar('Exit')           # this adds in a tool bar (actions linked to graphical icons)
        # self.toolbar.addAction(exitAction)               # this adds the exit action to the toolbar icon

        self.show()                        # this shows the GUI that's established



    def closeEvent(self, event):   # this is a function that displays a message box upon a quit event

        reply = QtGui.QMessageBox.question(self, 'Message', "Are you sure you want to quit?",
            QtGui.QMessageBox.Yes | QtGui.QMessageBox.No, QtGui.QMessageBox.No)

        if reply == QtGui.QMessageBox.Yes:
            event.accept()
        else:
            event.ignore()



def main():   # this is the "main" loop of the application, where the stuff described above is executed

    app = QtGui.QApplication(sys.argv)
    box = Example()
    sys.exit(app.exec_())

if __name__ == '__main__':   # this is the execution of the main loop function
    main()

