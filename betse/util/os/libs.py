#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level dynamically linked shared library facilities.

Caveats
----------
Shared library-specific logic is, by definition, operation system-specific logic
and hence poor form. Do so _only_ where necessary.
'''

# ....................{ IMPORTS                            }....................
from betse.util.type.types import type_check

# ....................{ GETTERS                            }....................
#FIXME: Add OS X support as well. Since OS X lacks the "ldd" command, doing
#so will require parsing the output of the equivalent OS X-specific "otool"
#command: e.g.,
#
#    # Under OS X.
#    $ otool -L h3dpost.x
#    h3dpost.x:
#       /opt/local/lib/libmpi.0.dylib (compatibility version 1.0.0, current version 1.0.0)
#       /usr/lib/libSystem.B.dylib (compatibility version 1.0.0, current version 111.1.3)

@type_check
def get_dependency_filenames(lib_filename: str) -> list:
    '''
    List of the absolute paths of all shared libraries dynamically linked to
    (and hence required at runtime) by the shared library with the passed path.

    Parameters
    ----------
    lib_filename : str
        Absolute or relative path of the shared library to inspect.

    Returns
    ----------
    list
        List of the absolute paths of all shared libraries required by this
        library.
    '''

    # Avoid circular import dependencies.
    from betse.util.os import oses
    from betse.util.type import regexes
    from betse.util.path import files
    from betse.util.path.command import pathables, runners

    #FIXME: Call die_unless_lib() instead, please.

    # If this library does *NOT* exist, raise an exception.
    files.die_unless_file(lib_filename)

    # If the current platform is Linux...
    if oses.is_linux():
        # String listing all libraries linked to by this library, captured from
        # the external "ldd" command.
        lib_ldd = runners.run_capturing_stdout(command_words=(
            writer_filename,
            '-help',
            'encoder=' + strs.shell_quote(codec_name),
        ))
        pass
        # regexes.iter_matches(
        #     text=, regex=, flags=regexes.FLAG_MULTILINE)

    #FIXME: Raise a human-readable exception.
    # Else, library inspection is currently unsupported on this platform.
    raise
