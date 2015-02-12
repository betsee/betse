#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level external process facilities.
'''

# ....................{ IMPORTS                            }....................
from betse.util.path import paths
import sys

# ....................{ GETTERS                            }....................
def get_current_basename() -> str:
    '''
    Get the basename of the executable originating the current process.
    '''
    # Since "sys.argv[0]" is either an absolute or relative path, get only such
    # path's basename.
    return paths.get_basename(sys.argv[0])

# --------------------( WASTELANDS                         )--------------------
