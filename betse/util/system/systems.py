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
import os, platform

# ....................{ TESTERS ~ os                       }....................
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

# ....................{ TESTERS ~ arch                     }....................
#FIXME: Shift to a new module "arches" in the same package.
#FIXME: Research us up. Is this *REALLY* the canonical means of determining
#this under Python3? (This seems fairly terrible, honestly. Is this actually
#cross-platform-portable?)

def is_32_bit():
    '''
    True if the current CPU architecture is 32-bit.
    '''
    return not is_64_bit()

def is_64_bit():
    '''
    True if the current CPU architecture is 64-bit.
    '''
    return 'PROCESSOR_ARCHITEW6432' in os.environ

# --------------------( WASTELANDS                         )--------------------
