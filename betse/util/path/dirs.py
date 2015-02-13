#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level directory facilities.

This module is named `dirs` rather than `dir` to avoid conflict with the `dir()`
builtin.
'''

# ....................{ IMPORTS                            }....................
from betse import metadata
from betse.util.exceptions import BetseExceptionDir
from betse.util.system import systems
from os import environ, path
import os

# ....................{ CONSTANTS                          }....................
#FIXME: Replace all existing calls to os.path.expanduser() with such constant.

HOME_DIR = path.expanduser('~')
'''
Absolute path of the current user's home directory.
'''

# ....................{ CONSTANTS ~ dot                    }....................
DOT_DIR = ''  # initialized below
'''
Absolute path of `betse`'s dot directory in the current user's home directory.

Specifically, under:

* Linux, this is `~/.betse/`. `betse` does *not* currently comply with the
  _XDG Base Directory Specification (e.g., `~/.local/share/betse`), which the
  principal authors of `betse` largely regard as unhelpful.
* OS X, this is `~/Library/Application Support/BETSE`.
* Windows, this is
  `C:\Documents and Settings\${User}\Application Data\BETSE`.

.. _XDG Base Directory Specification: http://standards.freedesktop.org/basedir-spec/basedir-spec-latest.html
'''

# If the current system is OS X, set such directory accordingly.
if systems.is_osx():
    DOT_DIR = path.join(
        HOME_DIR, 'Library', 'Application Support', metadata.SCRIPT_NAME_CLI)
# If the current system is Windows, set such directory accordingly.
elif systems.is_windows():
    DOT_DIR = path.join(environ['APPDATA'], metadata.NAME)
#FIXME: Explicitly assert POSIX compatibility here.
# Else, assume the current system is POSIX-compatible.
else:
    DOT_DIR = path.join(HOME_DIR, '.' + metadata.SCRIPT_NAME_CLI)

# ....................{ EXCEPTIONS                         }....................
def die_unless_found(dirname: str) -> None:
    '''
    Raise an exception unless the passed directory exists.
    '''
    assert isinstance(dirname, str), '"{}" not a string.'.format(dirname)
    if not is_dir(dirname):
        raise BetseExceptionDir(
            '"{}" not found, not a directory, or not readable.'.format(dirname))

def die_unless_parent_found(pathname: str) -> None:
    '''
    Raise an exception unless the parent directory of the passed path exists.
    '''
    assert isinstance(pathname, str), '"{}" not a string.'.format(pathname)
    die_unless_found(get_dirname(pathname))

# ....................{ TESTERS                            }....................
def is_dir(dirname: str) -> bool:
    '''
    True if the passed directory exists.
    '''
    assert isinstance(dirname, str), '"{}" not a string.'.format(dirname)
    return path.isdir(dirname)

# ....................{ GETTERS                            }....................
def get_dirname(pathname: str) -> str:
    '''
    Get the *dirname* (i.e., parent directory) of the passed path.
    '''
    assert isinstance(pathname, str), '"{}" not a string.'.format(pathname)
    return path.dirname(pathname)

# ....................{ JOINERS                            }....................
def join(*pathnames) -> str:
    '''
    Join the passed pathnames on the directory separator specific to the current
    operating system.

    This is a convenience function wrapping the standard `os.path.join()`
    function, provided to reduce the number of import statements required by
    other modules.
    '''
    return path.join(*pathnames)

# ....................{ MAKERS                             }....................
#FIXME: Replace all existing calls to os.makedirs() by calls to such function.

def make_unless_found(dirname: str) -> None:
    '''
    Create the passed directory if *not* found.

    All nonexistent parents of such directory will also be recursively created,
    mimicking the action of the conventional shell command `mkdir -p`.
    '''
    assert isinstance(dirname, str), '"{}" not a string.'.format(dirname)
    os.makedirs(dirname, exist_ok = True)

def make_parent_unless_found(*pathnames) -> None:
    '''
    Create the parent directory of each passed path if *not* found.
    '''
    for pathnames in pathnames:
        assert isinstance(pathname, str), '"{}" not a string.'.format(pathname)
        make_unless_found(get_dirname(pathname))

# --------------------( WASTELANDS                         )--------------------
# from betse.util.path import paths
