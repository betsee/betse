#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry
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
#FIXME: Implement me. The ideal implementation for each:
#
#1. Iteratively calls testers defined above for the three principal platforms
#   is_linux(), is_os_x(), and is_windows().
#2. If:
#   * is_os_x(), parse platform.mac_ver() output for the desired metadata.
#   * is_windows(), parse platform.win32_ver() output for the desired metadata.
#   * is_linux() and:
#     * If platform.linux_distribution() is defined, parse that function's
#       output for the desired metadata.
#     * Else, defer to:
#       * platform.release() for the version.
#       * platform.system() for the OS name.
#3. Else, defer to the same fallbacks as in the Linux case.

# def get_name() -> str:
#     '''
#     Human-readable name of the active Python interpreter's implementation (e.g.,
#     `CPython`, `PyPy`).
#     '''
#     return platform.python_implementation()
#
#
# def get_version() -> str:
#     '''
#     Human-readable `.`-delimited version specifier string of the current
#     platform.
#
#     Under:
#
#     * Linux, this is typically the version of the current Linux kernel (e.g.,
#       `4.1.15-gentoo-r1`).
#     * OS X, this is typically the version of the current Linux kernel (e.g.,
#       `4.1.15-gentoo-r1`).
#     '''
#     return platform.python_version()

# ....................{ GETTERS ~ metadata                 }....................
def get_metadata() -> OrderedDict:
    '''
    Get an ordered dictionary synopsizing the current platform.

    This function expands the metadata reported by the cross-platform function
    `platform.uname()` with additional diagnostics.
    '''

    # This dictionary.
    metadata = vars(platform.uname())

    #FIXME: Expand this dictionary here with additioral metadata. See:
    #    https://docs.python.org/3/library/platform.html
    #FIXME: Replace this dictionary's:
    #
    #* "system" key-value pair with a "name" key-value pair obtained by
    #  calling the get_name() function defined above.
    #* "release" key-value pair with a "version" key-value pair obtained by
    #  calling the get_version() function defined above.

    # Return this dictionary.
    return metadata
