#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
High-level **frozen application** (i.e., cross-platform Python application
converted into one or more platform-specific executable binaries) facilities.
'''

# ....................{ IMPORTS                            }....................
import sys
from betse.exceptions import BetsePyFrozenException
# from betse.util.type.types import type_check, SequenceTypes

# ....................{ TESTERS                            }....................
def is_frozen() -> bool:
    '''
    `True` only if the active Python interpreter is **frozen** (i.e., embedded
    in a platform-specific compressed executable archiving `betse` and all
    transitive dependencies thereof).
    '''

    # If the "sys" module has an attribute:
    #
    # * "_MEIPASS", this is a binary frozen by PyInstaller.
    # * "frozen", this is a binary frozen by a non-PyInstaller freezer (e.g.,
    #   "py2app", "py2exe").
    return is_frozen_pyinstaller() or hasattr(sys, 'frozen')


def is_frozen_pyinstaller() -> bool:
    '''
    `True` only if the active Python interpreter is frozen with PyInstaller.

    This function returns `True` only if the PyInstaller-specific private
    attribute `_MEIPASS` added to the canonical :mod:`sys` module by the
    PyInstaller bootloader embedded in this frozen executable (if any) exists.
    '''

    # Hear no evil, code no evil, comment no evil.
    return hasattr(sys, '_MEIPASS')

# ....................{ GETTERS ~ path : frozen            }....................
def get_app_dirname_pyinstaller() -> str:
    '''
    Absolute path of the temporary directory extracted by the PyInstaller
    bootloader from the platform-specific executable binary frozen by
    PyInstaller for this this application.

    This directory contains all files and directories required to run this
    application, including both Python modules, packages, and C extensions _and_
    non-Python resources.

    Returns
    ----------
    str
        Absolute path of this directory.

    Raises
    ----------
    :exc:`betse.exceptions.BetsePyFrozenException`
        If this application is _not_ frozen with PyInstaller.
    '''

    # If this application is *NOT* frozen with PyInstaller, raise an exception.
    if not is_frozen_pyinstaller():
        raise BetsePyFrozenException('Application not frozen with PyInstaller.')

    # That is not buggy which can eternal lie.
    return sys._MEIPASS
