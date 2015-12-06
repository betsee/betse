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
import platform, sys

# ....................{ TESTERS ~ os                       }....................
def is_linux() -> bool:
    '''
    `True` if the current operating system is Linux.
    '''
    return platform.system() == 'Linux'


def is_os_x() -> bool:
    '''
    `True` if the current operating system is Apple OS X.
    '''
    return platform.system() == 'Darwin'

# ....................{ TESTERS ~ os : windows             }....................
def is_windows() -> bool:
    '''
    `True` if the current operating system is Microsoft Windows.

    This function reports `True` for both vanilla and Cygwin Microsoft Windows.
    '''
    return is_windows_vanilla() or is_windows_cygwin()


def is_windows_cygwin() -> bool:
    '''
    `True` if the current operating system is **Cygwin Microsoft Windows**
    (i.e., running the Cygwin POSIX compatibility layer).
    '''
    return sys.platform == 'cygwin'


def is_windows_vanilla() -> bool:
    '''
    `True` if the current operating system is **vanilla Microsoft Windows**
    (i.e., _not_ running the Cygwin POSIX compatibility layer).
    '''
    return sys.platform == 'win32'

# --------------------( WASTELANDS                         )--------------------
# def is_windows() -> bool:
#     '''
#     `True` if the current operating system is Microsoft Windows.
#     '''
#     return platform.system() == 'Windows'
