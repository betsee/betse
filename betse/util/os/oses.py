#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
High-level operating system (OS) facilities.

Caveats
----------
Operating system-specific logic is poor form and should be leveraged only where
necessary.
'''

# ....................{ IMPORTS                            }....................
import platform, sys
from collections import OrderedDict

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

# ....................{ GETTERS                            }....................
def get_metadata() -> OrderedDict:
    '''
    Get an ordered dictionary synopsizing the current system.

    This function expands the metadata reported by the cross-platform function
    `platform.uname()` with additional diagnostics.
    '''
    # This dictionary.
    metadata = vars(platform.uname())

    #FIXME: Expand this dictionary here with additioral metadata. See:
    #    https://docs.python.org/3/library/platform.html

    # Get such dictionary.
    return metadata

# --------------------( WASTELANDS                         )--------------------
# def is_windows() -> bool:
#     '''
#     `True` if the current operating system is Microsoft Windows.
#     '''
#     return platform.system() == 'Windows'
