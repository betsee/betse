#!/usr/bin/env python3

import sys
from PySide import QtGui
from PySide import QtCore
import floatvalidator

# sys.path.append('/home/pietakio/BioE/test')

class DoubleEdit(QtGui.QLineEdit):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setValidator(floatvalidator.float_validator)

class MainApp(QtGui.QMainWindow):

    def __init__(self):
        super().__init__()
        self.initUi()

    def initUi(self):
        self.mainframe = QtGui.QWidget()
        self.lab_1 = QtGui.QLabel(self.tr('Na+ permeability'))
        self.textbox_1 = DoubleEdit('10e-9')
        self.OKbtn = QtGui.QPushButton('Save Values')

        self.OKbtn.clicked.connect(self.setparams)

        self.lab_2 = QtGui.QLabel(self.tr('K+ permeability'))
        self.textbox_2 = DoubleEdit('10e-9')
        vboxlay = QtGui.QVBoxLayout()
        vboxlay.addWidget(self.lab_1)
        vboxlay.addWidget(self.textbox_1)
        vboxlay.addWidget(self.lab_2)
        vboxlay.addWidget(self.textbox_2)
        vboxlay.addWidget(self.OKbtn)
        self.mainframe.setLayout(vboxlay)
        self.setCentralWidget(self.mainframe)

    def setparams(self):
        try:
            P_Na = float(self.textbox_1.text())
            K_Na = float(self.textbox_2.text())
            print(P_Na, K_Na)

        except:
            print("format is invalid")
            self.textbox_1.setText('10e-9')
            self.textbox_2.setText('10e-9')


def main():   # this is the "main" loop of the application, where the stuff described above is executed
    app = QtGui.QApplication(sys.argv)
    appview = MainApp()
    appview.setWindowTitle('YoLo Babe!')
    appview.show()
    sys.exit(app.exec_())

if __name__ == '__main__':   # this is the execution of the main loop function
    main()