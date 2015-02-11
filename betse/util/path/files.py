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
from betse.util.path import dirs
import os

# ....................{ CONSTANTS                          }....................
#FIXME: Actually use this. More work than we currently care to invest.

DEFAULT_CONFIG_FILE = dirs.join(dirs.DOT_DIR, 'config.yaml')
'''
Absolute path of the default user-specific file with which `betse` configures
application-wide behaviour (e.g., log settings).
'''

DEFAULT_LOG_FILE = dirs.join(dirs.DOT_DIR, 'debug.log')
'''
Absolute path of the default user-specific file to which `betse` logs messages.
'''

# ....................{ REMOVERS                           }....................
def remove(filename: str) -> None:
    '''
    Remove the passed file.

    If such file does *not* exist, an exception is raised.
    '''
    assert isinstance(filename, str), '"{}" not a string.'.format(filename)
    os.remove(filename)

# --------------------( WASTELANDS                         )--------------------
