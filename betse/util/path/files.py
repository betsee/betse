#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level non-directory file facilities.

This module is named `files` rather than `file` to avoid conflict with the stock
`files` class.
'''

# ....................{ IMPORTS                            }....................
import os

# ....................{ REMOVERS                           }....................
def remove(filename: str) -> None:
    '''
    Remove the passed file.

    If such file does *not* exist, an exception is raised.
    '''
    assert isinstance(filename, str), '"{}" not a string.'.format(filename)
    os.remove(filename)

# --------------------( WASTELANDS                         )--------------------
