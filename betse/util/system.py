#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level operating system facilities.

Caveats
----------
Operating system-specific functions (e.g., `is_windows()`) are considered poor
form. Call such functions only when absolutely necessary.
'''

# ....................{ IMPORTS                            }....................
import platform

# ....................{ OUTPUTTERS                         }....................
def is_linux() -> bool:
    '''
    True if the current operating system is Linux.
    '''
    return platform.system() == 'Linux'

def is_osx() -> bool:
    '''
    True if the current operating system is Apple OS X.
    '''
    return platform.system() == 'Darwin'

def is_windows() -> bool:
    '''
    True if the current operating system is Microsoft Windows.
    '''
    return platform.system() == 'Windows'

# --------------------( WASTELANDS                         )--------------------
