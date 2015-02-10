#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level external process facilities.
'''

# ....................{ IMPORTS                            }....................
import sys

# ....................{ GETTERS                            }....................
def get_current_basename() -> str:
    '''
    Get the basename of the executable originating the current process.
    '''
    return sys.argv[0]

# --------------------( WASTELANDS                         )--------------------
