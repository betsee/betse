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
from betse.util import system
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
if system.is_osx():
    DOT_DIR = path.join(
        HOME_DIR, 'Library', 'Application Support', metadata.SCRIPT_NAME_CLI)
# If the current system is Windows, set such directory accordingly.
elif system.is_windows():
    DOT_DIR = path.join(environ['APPDATA'], metadata.NAME)
#FIXME: Explicitly assert POSIX compatibility here.
# Else, assume the current system is POSIX-compatible.
else:
    DOT_DIR = path.join(HOME_DIR, '.' + metadata.SCRIPT_NAME_CLI)

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

# --------------------( WASTELANDS                         )--------------------
