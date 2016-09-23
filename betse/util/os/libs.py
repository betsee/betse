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
from betse.exceptions import BetseLibException, BetseOSException
from betse.util.type.types import type_check, GeneratorType

# ....................{ CONSTANTS                          }....................
KERNEL_NAME_TO_LIB_FILETYPES = {
    # *ALL* Linux distributions.
    'Linux':   {'so',},

    # OS X and iOS. This currently only includes Mach-O shared libraries.
    # Macho-O bundles (i.e., loadable modules) are interchangeable with Macho-O
    # shared libraries for some but *NOT* all operations and hence omitted.
    'Darwin':  {'dylib',},

    # Windows. Technically, other Windows-specific shared library formats exist
    # (e.g., "ocx" for libraries containing ActiveX controls and "drv" for
    # legacy system drivers). Since it's unclear whether anyone still cares
    # about these formats, they are currently omitted.
    'Windows': {'dll',},
}
'''
Dictionary mapping from the name of each possible platform returned by the
:func:`betse.util.os.kernels.get_name` function to the set of all undotted
filetypes signifying shared libraries supported by this platform.
'''

# ....................{ EXCEPTIONS                         }....................
def die_unless_lib(pathname: str) -> None:
    '''
    Raise an exception unless the passed path is that of an existing shared
    library supported by the current platform.

    Parameters
    ----------
    lib_filename : str
        Absolute or relative path of the shared library to inspect.
    '''

    if not is_lib(pathname):
        raise BetseLibException(
            '"{}" not a shared library.'.format(pathname))

# ....................{ TESTERS                            }....................
@type_check
def is_lib(pathname: str) -> bool:
    '''
    `True` only if the passed path is that of an existing shared library
    specific to the current platform (e.g., suffixed by `.so` on Linux).

    This function returns `False` if either:

    * This path is _not_ that of an existing file.
    * This pathname has no filetype.
    * This pathname's filetype is _not_ that of a shared library supported by
      the current platform.

    Parameters
    ----------
    pathname : str
        Absolute or relative path of the shared library to inspect.

    Returns
    ----------
    bool
        `True` only if this path is that of an existing shared library.
    '''

    # Avoid circular import dependencies.
    from betse.util.path import files, paths
    from betse.util.os import kernels

    # Filetype of this pathname if any or None otherwise.
    filetype = paths.get_filetype_undotted_or_none(pathname)

    # Name of this platform's low-level kernel (e.g., "Linux", "Darwin").
    kernel_name = kernels.get_name()

    # This path is that of an existing shared library if and only if...
    return (
        # This pathname has a filetype.
        filetype is not None and
        # This filetype is that of a shared library for this platform.
        filetype in KERNEL_NAME_TO_LIB_FILETYPES[kernel_name] and
        # This path is that of an existing file. Since testing this requires a
        # filesystem lookup, this is the slowest test and thus deferred.
        files.is_file(pathname)
    )

# ....................{ GENERATORS                         }....................
#FIXME: Add OS X support as well. Since OS X lacks the "ldd" command, doing
#so will require parsing the output of the equivalent OS X-specific "otool"
#command: e.g.,
#
#    # Under OS X.
#    $ otool -L h3dpost.x
#    h3dpost.x:
#       /opt/local/lib/libmpi.0.dylib (compatibility version 1.0.0, current version 1.0.0)
#       /usr/lib/libSystem.B.dylib (compatibility version 1.0.0, current version 111.1.3)

def iter_linked_lib_filenames(lib_filename: str) -> GeneratorType:
    '''
    Generator iteratively yielding the 2-tuple of the basename and absolute path
    of each shared library dynamically linked to (and hence required at runtime
    by) the shared library with the passed path.

    Parameters
    ----------
    lib_filename : str
        Absolute or relative path of the shared library to inspect.

    Returns
    ----------
    GeneratorType
        Generator iteratively yielding the 2-tuple
        `(linked_lib_basename, linked_lib_pathname`), where:
        * `linked_lib_basename` is the basename of a shared library dynamically
          linked to the shared library with the passed path.
        * `linked_lib_pathname` is the absolute pathname of the same library.
    '''

    # Avoid circular import dependencies.
    from betse.util.os import oses
    from betse.util.path.command import runners
    from betse.util.type import regexes, strs

    # If this library does *NOT* exist, raise an exception.
    die_unless_lib(lib_filename)

    # If the current platform is Linux...
    if oses.is_linux():
        # Tuple of all shell words comprising the "ldd" command to be run.
        ldd_command_words = ('ldd', strs.shell_quote(lib_filename))

        # String listing all libraries linked to by this library captured from
        # stdout printed by this command. For example, when passed the absolute
        # path of the file defining the "numpy.core.multiarray" C extension,
        # this command prints stdout resembling:
        #
        # multiarray.cpython-34.so:
        # 	linux-vdso.so.1 (0x00007fff049ca000)
        # 	libopenblas_threads.so.0 => /usr/lib64/libopenblas_threads.so.0 (0x00007f9af3ad5000)
        # 	libm.so.6 => /lib64/libm.so.6 (0x00007f9af37d9000)
        # 	libpython3.4.so.1.0 => /usr/lib64/libpython3.4.so.1.0 (0x00007f9af3357000)
        # 	libc.so.6 => /lib64/libc.so.6 (0x00007f9af2fc1000)
        # 	/lib64/ld-linux-x86-64.so.2 (0x000055e2dd7f3000)
        # 	libpthread.so.0 => /lib64/libpthread.so.0 (0x00007f9af2da5000)
        # 	libdl.so.2 => /lib64/libdl.so.2 (0x00007f9af2ba0000)
        # 	libutil.so.1 => /lib64/libutil.so.1 (0x00007f9af299d000)
        ldd_stdout = runners.run_capturing_stdout(ldd_command_words)

        # For each line containing a "=>"-delimited basename and absolute path
        # pair and hence ignoring both pseudo-libraries that do *NOT* actually
        # exist (e.g., "linux-vdso") or ubiquitous libraries required by the
        # dynamic linking mechanism and hence guaranteed to *ALWAYS* exist
        # (e.g., "ld-linux-x86-64")...
        for line_match in regexes.iter_matches_line(
            text=ldd_stdout,
            regex=r'^\s+(\S+)\s+=>\s+(\S+)\s+\(0x[0-9a-fA-F]+\)$',
        ):
            # Yield the 2-tuple corresponding exactly to the match groups
            # captured by this match.
            yield line_match.groups()

    # Else, library inspection is currently unsupported on this platform.
    raise BetseOSException(
        'Shared library inspection currently unsupported under {}.'.format(
            oses.get_name()))
