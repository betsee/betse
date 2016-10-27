#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level **path archival** (i.e., compression and decompression of both
directories _and_ non-directory files) facilities.
'''

#FIXME: Actually leverage this submodule to transparently archive pickles.

# ....................{ IMPORTS                            }....................
from betse.exceptions import BetseArchiveException
from betse.util.type.types import type_check
from io import BufferedIOBase

# Unlike most archive-specific stdlib modules, this module is unconditionally
# installed by *ALL* Python installations and hence safely importable here.
from zipfile import ZipFile, ZIP_DEFLATED

# ....................{ CONSTANTS ~ public                 }....................
# Initialized below by the init() function.
ARCHIVE_FILETYPES = set()
'''
Set of all archive filetypes supported by this submodule.
'''

# ....................{ EXCEPTIONS                         }....................
@type_check
def die_unless_filetype_archive(pathname: str) -> None:
    '''
    Raise an exception unless the passed pathname is suffixed by a filetype
    corresponding to a supported archive format.

    Parameters
    ----------
    pathname : str
        Relative or absolute path of the file to be validated.

    See Also
    ----------
    :func:`is_filetype_archive`
        Further details.
    '''

    # If this pathname is *NOT* suffixed by a supported archive filetype...
    if not is_filetype_archive(pathname):
        # Avoid circular import dependencies.
        from betse.util.path import paths
        from betse.util.py import modules
        from betse.util.type import strs

        # Filetype of this pathname if any or raise an exception otherwise.
        filetype = paths.get_filetype_undotted(pathname)

        # Fully-qualified name of the optional stdlib module required to write
        # archives of this filetype if this is a recognized (but unsupported)
        # archive filetype or None otherwise.
        module_name = _ARCHIVE_FILETYPE_TO_MODULE_NAME.get(filetype, None)

        # If this is a recognized archive filetype, this module *CANNOT* have
        # been installed with Pythoni, implying the corresponding third-party
        # library headers to have been unavailable at Python installation time.
        if module_name is not None:
            # Exception message to be raised describing this constraint.
            exception_message = (
                'Module "{}" not found '
                '(typically due to third-party library headers '
                'required by this module being unavailable at '
                'Python installation time).'.format(module_name))

            # Raise the expected exception.
            modules.die_unless_module(module_name, exception_message)

            # The above call should have raised an exception. For safety,
            # unconditionally raise another exception as a sentinel.
            raise ValueError('That which should not be has become.')
        #  Else, this is an unrecognized archive filetype.
        else:
            # Human-readable string listing all supported archive filetypes.
            archive_filetypes = strs.join_as_disconjunction_double_quoted(
                *ARCHIVE_FILETYPES)

            # Raise an exception embedding this string.
            raise BetseArchiveException(
                'Filename "{}" not suffixed by a '
                'supported archive format (i.e., {}).'.format(
                    pathname, archive_filetypes))

# ....................{ TESTERS                            }....................
@type_check
def is_filetype_archive(pathname: str) -> bool:
    '''
    `True` only if the passed pathname is suffixed by a filetype corresponding
    to a supported archive format.

    `True` is returned only if the **rightmost filetype** (i.e., substring
    following the last `.` character) of this pathname is in the
    :data:`ARCHIVE_FILETYPES` set global. All preceding filetypes of this
    pathname are ignored.

    Parameters
    ----------
    pathname : str
        Relative or absolute path of the file to be tested.

    Returns
    ----------
    bool
        `True` only if this pathname is suffixed by an archive filetype.
    '''

    # Avoid circular import dependencies.
    from betse.util.path import paths

    # True only if this pathname is suffixed by an archive filetype.
    return paths.is_filetype_undotted_in(
        pathname=pathname, filetypes_undotted=ARCHIVE_FILETYPES)

# ....................{ WRITERS                            }....................
#FIXME: Unit test us up.
@type_check
def write_bytes(filename: str) -> BufferedIOBase:
    '''
    Open and return a filehandle suitable for writing the binary archive file
    with the passed filename.

    This function returns a :class:`file`-like object suitable for use wherever
    the :func:`open` builtin is callable (e.g., in `with` statements).

    Parameters
    ----------
    filename : str
        Relative or absolute path of the binary archive file to be written,
        whose **rightmost filetype** (i.e., substring suffixing the last `.`
        character in this pathname) _must_ be that of a supported archive
        format.

    Returns
    ----------
    BufferedIOBase
        `file`-like object encapsulating this opened file.

    Raises
    ----------
    BetseArchiveException
        If this filename is _not_ suffixed by a supported archive filetype.
    BetsePathException
        If this filename is that of an existing file or directory.
    '''

    # Avoid circular import dependencies.
    from betse.util.path import dirs, paths

    # Raise an exception unless this filename has a supported archive filetype
    # *AND* is not already an existing file or directory.
    die_unless_filetype_archive(filename)
    paths.die_if_path(filename)

    # Create the parent directory of this file if needed.
    dirs.make_parent_unless_dir(filename)

    # Filetype of this pathname. By the above validation, this filetype is
    # guaranteed to both exist and be a supported archive filetype. Hence, no
    # additional validation is required.
    filetype = paths.get_filetype_undotted(filename)

    # Callable writing archives of this filetype.
    writer = _ARCHIVE_FILETYPE_TO_WRITER[filetype]

    # Open and return a filehandle suitable for writing this archive.
    return writer(filename)

# ....................{ WRITERS ~ private                  }....................
def _write_bytes_bz2(filename: str) -> BufferedIOBase:
    '''
    Open and return a filehandle suitable for writing the binary bzip-archived
    file with the passed filename.

    This private function is intended to be called _only_ by the public
    :func:`write_bytes` function.
    '''

    # This optional stdlib module is guaranteed to exist and hence be safely
    # importable here, due to the above die_unless_filetype_archive() call.
    from bz2 import BZ2File

    # Open and return a filehandle suitable for e(x)clusively writing this file.
    return BZ2File(filename, mode='xb')


def _write_bytes_gz(filename: str) -> BufferedIOBase:
    '''
    Open and return a filehandle suitable for writing the binary gzip-archived
    file with the passed filename.

    This private function is intended to be called _only_ by the public
    :func:`write_bytes` function.
    '''

    # This optional stdlib module is guaranteed to exist and hence be safely
    # importable here, due to the above die_unless_filetype_archive() call.
    from gzip import GzipFile

    # Open and return a filehandle suitable for e(x)clusively writing this file.
    return GzipFile(filename, mode='xb')


def _write_bytes_xz(filename: str) -> BufferedIOBase:
    '''
    Open and return a filehandle suitable for writing the binary LZMA-archived
    file with the passed filename.

    This private function is intended to be called _only_ by the public
    :func:`write_bytes` function.
    '''

    # This optional stdlib module is guaranteed to exist and hence be safely
    # importable here, due to the above die_unless_filetype_archive() call.
    from lzma import LZMAFile

    # Open and return a filehandle suitable for e(x)clusively writing this file.
    return LZMAFile(filename, mode='xb')


def _write_bytes_zip(filename: str) -> BufferedIOBase:
    '''
    Open and return a filehandle suitable for writing the binary zip-archived
    file with the passed filename.

    This private function is intended to be called _only_ by the public
    :func:`write_bytes` function.
    '''

    # Open and return a filehandle suitable for e(x)clusively writing this file.
    return ZipFile(
        filename=filename,
        mode='xb',

        # Compress this file with the optional "zlib" stdlib module, which is
        # guaranteed to exist and hence be safely indirectly importable here due
        # to the above die_unless_filetype_archive() call. By default,
        # compression defaults to the "ZIP_STORED" constant disabling
        # compression -- which is frankly insane, but ultimately unsurprising.
        compression=ZIP_DEFLATED,
    )

# ....................{ CONSTANTS ~ private                }....................
_ARCHIVE_FILETYPE_TO_MODULE_NAME = {
    'bz2': 'bz2',
    'gz':  'gzip',
    'xz':  'lzma',

    # The low-level optional "zlib" module provides the compression
    # functionality required by the high-level mandatory "zipfile" module. While
    # the latter is technically usable without the former, the resulting
    # archives will be uncompressed and hence effectively useless.
    'zip': 'zlib',
}
'''
Dictionary mapping from each archive filetype supported by this submodule to
the fully-qualified name of the optional stdlib module required to write
archives of this filetype.
'''


_ARCHIVE_FILETYPE_TO_WRITER = {
    'bz2': _write_bytes_bz2,
    'gz':  _write_bytes_gz,
    'xz':  _write_bytes_xz,
    'zip': _write_bytes_zip,
}
'''
Dictionary mapping from each archive filetype supported by this submodule to
the private function of this submodule writing archives of this filetype.
'''

# ....................{ INITIALIZERS                       }....................
def init() -> None:
    '''
    Initialize this submodule.

    This function defines all currently undefined global submodule variables.
    '''

    # Avoid circular import dependencies.
    from betse.util.py import modules

    # Globals initialized below.
    global ARCHIVE_FILETYPES

    # For each supported archive filetype and the fully-qualified name of the
    # optional stdlib module required to write archives of this filetype....
    for filetype, module_name in _ARCHIVE_FILETYPE_TO_MODULE_NAME.items():
        # If this module is installed...
        if modules.is_module(module_name):
            # Add this filetype to the set of all supported archive filetypes.
            ARCHIVE_FILETYPES.add(filetype)

# ....................{ MAIN                               }....................
# Initialize this submodule.
init()
