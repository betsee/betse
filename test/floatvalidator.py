#!/usr/bin/env python3
from PySide import QtGui
import re

class FloatValidator(QtGui.QValidator):
    '''
    This is...
    '''

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._float_re = re.compile(r'(([+-]?\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?)')

    def is_valid_float_string(self,string):
        match = self._float_re.search(string)
        return match.groups()[0] == string if match else False
        # self._float_re.match(string) is None

    def validate(self, string, position):
        if self.is_valid_float_string(string):
            return self.State.Acceptable
        elif string == "" or string[position-1] in 'e.-+':   #FIXME this accepts eee +++ --- as valid!
            return self.State.Intermediate
        return self.State.Invalid

    def fixup(self, text):
        match = self._float_re.search(text)
        return match.groups()[0] if match else ""


float_validator = FloatValidator()
