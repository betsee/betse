#!/usr/bin/env python3

import sys
from PySide import QtCore
from PySide import QtGui



class Example(QtGui.QWidget):

    def __init__(self):
        super().__init__()

        self.initUi()


    def initUi(self):

        self.setWindowTitle("Yolo")
        self.setGeometry(300,300,500,500)
        self.centre()
        icon = QtGui.QIcon()
        icon.addFile("/home/pietakio/BioE/test/android.gif")
        #self.setWindowIcon(icon)
        button = QtGui.QPushButton("Push", self)
        button.setIcon(icon)
        button.setToolTip("This is a robo button")
        button.move(250,250)
        qbtn = QtGui.QPushButton("Quit",self)
        qbtn.setIcon(QtGui.QIcon("plant.png"))
        qbtn.clicked.connect(QtCore.QCoreApplication.instance().quit)
        self.show()

    def centre(self):

        qr = self.frameGeometry()
        cp = QtGui.QDesktopWidget().availableGeometry().center()
        qr.moveCenter(cp)
        self.move(qr.topLeft())

    def closeEvent(self, event):

        reply = QtGui.QMessageBox.question(self, 'Message', "Are you sure you want to quit?",
            QtGui.QMessageBox.Yes | QtGui.QMessageBox.No, QtGui.QMessageBox.No)

        if reply == QtGui.QMessageBox.Yes:
            event.accept()
        else:
            event.ignore()



def main():

    app = QtGui.QApplication(sys.argv)
    box = Example()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()

