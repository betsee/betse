#!/usr/bin/env python3
# GUI training

import sys
from PySide.QtCore import *
from PySide.QtGui import *
import math
import time

# Create a Qt application
#app = QApplication(sys.argv)
# Create a Label and show it
#label = QLabel("<font color=blue size=40>Love Works</font>")
#label.show()
# Enter Qt application main loop
#app.exec_()
#sys.exit()

# Another one

#class Form(QDialog):
#    def __init__(self, parent=None):   # FIXME maybe Sess can help explain these...
#        super().__init__(parent)
#        self.setWindowTitle("Love Bot Buttons")
#        self.edit = QLineEdit("Write name here")     # create a text input box widget
#        self.button = QPushButton("Show Greeting")   # create a button widget
#        layout = QVBoxLayout()                       # call an instance of the QVBoxLayout class
#        layout.addWidget(self.edit)                  # add the widgets to the layout
#        layout.addWidget(self.button)
#        self.setLayout(layout)                       # establish the layout as *the* layout we just made
#        self.button.clicked.connect(self.greetings)

#    def greetings(self):
 #       print("Hello %s" % self.edit.text())

#if __name__ == '__main__':   # FIXME I have no clue what that does!
#app = QApplication(sys.argv)   # create the application
#form = Form()
#form.show()
   # sys.exit(app.exec_())  # FIXME I don't get why this is inside the sys.exit command???

#sys.exit()
#app.exec_()
#sys.exit

# Now the calculator

# class Form(QDialog):
#
#      def __init__(self,parent=None):
#         super(Form,self).__init__(parent)
#
#         self.browser = QTextBrowser()
#         self.lineedit = QLineEdit("Type an expression & press enter")
#         self.lineedit.selectAll()
#
#         layout = QVBoxLayout()
#         layout.addWidget(self.browser)
#         layout.addWidget(self.lineedit)
#         self.setLayout(layout)
#
#         self.lineedit.setFocus()   # means this will be where the cursor starts
#         #self.connect(self.lineedit,SIGNAL("returnPressed()"),self.updateUi)
#         self.lineedit.returnPressed.connect(self.updateUi)
#         self.setWindowTitle("Calculate")
#
#     def updateUi(self):
#         try:
#             text=self.lineedit.text()
#             #self.browser.append("%s" % (text, math.eval(text)))
#             self.browser.append(text)
#
#         except:
#             self.browser.append("is invalid")
#
# app = QApplication(sys.argv)
# form = Form()
# form.show()
# app.exec_()







