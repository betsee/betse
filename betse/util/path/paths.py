#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level path facilities specific to neither directories nor non-directory
files.

This module is named `paths` rather than `path` to avoid conflict with the stock
`path` module of the `os` package.
'''

# ....................{ IMPORTS                            }....................
from os import path

# ....................{ GETTERS                            }....................
def get_basename(pathname: str) -> str:
    '''
    Get the *basename* (i.e., last component) of the passed path.
    '''
    assert isinstance(pathname, str), '"{}" not a string.'.format(pathname)
    return path.basename(pathname)

# ....................{ MAKERS                             }....................
def init() -> None:
    '''
    Validate core directories and files required at program startup.

    This function automatically creates non-existent paths where feasible and
    otherwise raises exceptions on such paths *not* being found or *not* having
    correct metadata (e.g., permissions).

    Such paths are required by both the CLI and GUI interfaces for `betse`. To
    support caller-specific exception handling, this function *must* be manually
    called sufficiently early in such startup.
    '''
    # Defer such imports to avoid circular dependencies.
    from betse.util.path import dirs

    #FIXME: Uncomment this after we ascertain the correct path for DATA_DIRNAME,
    #which we should probably get around to!

    # Make betse's top-level dot directory if not found.
    # dirs.die_unless_found(dirs.DATA_DIRNAME)

    # Make betse's top-level dot directory if not found.
    dirs.make_unless_found(dirs.DOT_DIRNAME)

# ....................{ REMOVERS                           }....................
# def remove(filename: str) -> None:
#     '''
#     Remove the passed file.
#
#     If such file does *not* exist, an exception is raised.
#     '''
#     assert isinstance(filename, str), '"{}" not a string.'.format(filename)
#     os.remove(filename)

# --------------------( WASTELANDS                         )--------------------
