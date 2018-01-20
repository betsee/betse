#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level **path archival** (i.e., compression and decompression of both
directories _and_ non-directory files) facilities.
'''

#FIXME: Add a new convert_format_to_filetype() function defined as follows
#
#    @type_check
#    def convert_format_to_filetype(format: str) -> StrOrNoneType:
#
#        #FIXME: Implement this format. Basically, we need to define an additional
#        #private sequence global -- say:
#        #
#        #    _ARCHIVE_FILETYPES_PREFERRED = ('bosc', 'lzo', 'xz', 'bz2', 'gz')
#        #
#        #Given that, we then need to perform the intersection of that sequence with
#        #the "ARCHIVE_FILETYPES" set, producing a new sequence containing only the
#        #elements of that sequence supported by the current system. Lastly, we then
#        #need to return the *FIRST* element of that sequence as the implementation
#        #of this format. For safety, note that this intersection could
#        #technically be empty, in which case None should be silently returned.
#        if format == 'auto':
#            ...  # something something
#        elif format == 'none':
#            return None
#        elif format in ARCHIVE_FILETYPES:
#            return format
#        #FIXME: It'd be great to raise a more comprehensive exception, but
#        #this would certainly do for the moment.
#        else:
#            raise BetseArchiveException(
#                'Archive format "{}" unsupported.'.format(format))
#FIXME: Refactor the "FILE HANDLING" section of our YAML file to resemble:
#
#path:
#  seed:
#    object:
#      file: SEEDS/world_1.betse  # File containing all cell cluster objects.
#
#  init:             # Paths for initialization saving and loading.
#    object:
#      file: INITS/init_1.betse   # File containing all initialization objects.
#    export:
#      directory: RESULTS/init_1  # Directory containing all initialization
#                                 # exports (e.g., plots, animations).
#
#  sim:
#    object:
#      file: SIMS/sim_1.betse    # filename to load/save sim run
#    export:
#      directory: RESULTS/sim_1  # directory to load/save data, plots, and animations for sim
#
#  object:
#    compression: auto  # Supported object compression formats include:
#                       # "auto" (smart compression), "none" (no compression),
#                       # "gz" (fast minimal compression),
#                       # "bz2" (medium compression), and
#                       # "xz" (slow maximal compression).
#FIXME: Pass the current value of the above "path/object/compression" key to the
#convert_format_to_filetype() function, producing a filetype suitable for use in
#pickling and unpickling compressed archives of that compression type.

#FIXME: Add support for additional non-stdlib archive packages, including:
#
#* https://github.com/Blosc/python-blosc, "...a Python wrapper for the extremely
#  fast Blosc compression library."
#* https://github.com/Iotic-Labs/py-lz4framed and
#  https://github.com/python-lz4/python-lz4, competing Python bindings for the
#  LZ4 compression library.
#* https://github.com/sergey-dryabzhinsky/python-zstd and
#  https://github.com/indygreg/python-zstandard, competing Python bindings for
#  the ZSTD compression library.
#
#While all of the above packages are released under BSD-compatible licenses,
#note that the following non-stdlib archive packages are strictly constrained by
#GPL licenses and *MUST* thus be avoided:
#
#* https://github.com/jd-boyd/python-lzo.
#
#Since the already supported XZ (LZMA) format arguably satisfies the ideal
#tradeoff between compression speed and ratio for our purposes, there may not
#exist much incentive to actually implement support for additional non-stdlib
#formats. Nonetheless, the future exists and it is always coming.

# ....................{ IMPORTS                            }....................
from io import BufferedIOBase

from betse.exceptions import BetseArchiveException
from betse.util.type.types import type_check

# ....................{ CONSTANTS ~ public                 }....................
# Initialized below by the init() function.
ARCHIVE_FILETYPES = set()
'''
Set of all archive filetypes supported by this submodule.

Tradeoffs
----------
Most compression algorithms succumb to the typical tradeoff between time and
space complexity. Specifically, as the **(de)compression speed** (i.e., duration
of time required to compress or decompress a given sequence of bytes) of a
compression algorithm increases, the **compression ratio** (i.e., filesize of
the output compressed archive to the size in bytes of the uncompressed input)
of that algorithm typically increases. The principal exception to this heuristic
is BZIP2, which runs slower and compresses less than competing algorithms.

Most benchmarks qualitatively profile the algorithms supported by this submodule
as follows (_in increasing order of compression ratio_):

* GZIP, exhibiting the worst compression ratio and fastest compression speed.
* BZIP2, exhibiting average compression ratio and slowest compression speed.
* XZ (effectively equivalent to LZMA), exhibiting the best compression ratio and
  average compression speed.

For general-purpose use, the **XZ** algorithm arguably strikes ideal balance
between compression ratio and speed. Most algorithms exhibit similar
decompression speed, which is thus ignorable for comparison purposes.

Caveats
----------
This set contains neither the ``tar`` nor ``zip`` container filetypes, which
remain intentionally unsupported by this submodule. Neither the standard
:mod:`tarfile` module supporting ``tar`` files nor the standard :mod:`zipfile`
module supporting ``zip`` files provide a :class:`file`-like API for reading and
writing arbitrary bytes or characters to and from these files, the primary use
case and hence mandate of this submodule.

While wrapping the container-oriented APIs provided by these modules with
:class:`file`-like APIs is feasible, there exists no demonstrable incentive to
do so. All compression formats supported by these modules are conveniently
already supported by :class:`file`-like APIs implemented by other standard
modules (e.g., the :class:`bz2.BZ2File` and :class:`gzip.GzipFile` classes).
'''

# ....................{ EXCEPTIONS                         }....................
@type_check
def die_unless_filetype(pathname: str) -> None:
    '''
    Raise an exception unless the passed pathname is suffixed by a filetype
    corresponding to a supported archive format.

    Parameters
    ----------
    pathname : str
        Relative or absolute path of the file to be validated.

    See Also
    ----------
    :func:`is_filetype`
        Further details.
    '''

    # If this pathname is *NOT* suffixed by a supported archive filetype...
    if not is_filetype(pathname):
        # Avoid circular import dependencies.
        from betse.util.path import pathnames
        from betse.util.type import modules
        from betse.util.type.text import strs

        # Filetype of this pathname if any or raise an exception otherwise.
        filetype = pathnames.get_filetype_undotted(pathname)

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
def is_filetype(pathname: str) -> bool:
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
    from betse.util.path import pathnames

    # True only if this pathname is suffixed by an archive filetype.
    return pathnames.is_filetype_undotted_in(
        pathname=pathname, filetypes_undotted=ARCHIVE_FILETYPES)

# ....................{ IO                                 }....................
@type_check
def read_bytes(filename: str) -> BufferedIOBase:
    '''
    Open and return a filehandle suitable for reading the binary archive file
    with the passed filename.

    This function returns a :class:`file`-like object suitable for use wherever
    the :func:`open` builtin is callable (e.g., in ``with`` statements).

    Parameters
    ----------
    filename : str
        Relative or absolute path of the binary archive file to be read, whose
        **rightmost filetype** (i.e., substring suffixing the last ``.``
        character in this pathname) *must* be that of a supported archive
        format.

    Returns
    ----------
    BufferedIOBase
        :class:`file`-like object encapsulating this opened file.

    Raises
    ----------
    BetseArchiveException
        If this filename is *not* suffixed by a supported archive filetype.
    BetseFileException
        If this filename is *not* that of an existing file.
    '''

    # Avoid circular import dependencies.
    from betse.util.path import files, pathnames

    # Raise an exception unless this filename has a supported archive filetype
    # *AND* is an existing file.
    die_unless_filetype(filename)
    files.die_unless_file(filename)

    # Filetype of this pathname. By the above validation, this filetype is
    # guaranteed to both exist and be a supported archive filetype. Hence, no
    # additional validation is required.
    filetype = pathnames.get_filetype_undotted(filename)

    # Callable reading archives of this filetype.
    reader = _ARCHIVE_FILETYPE_TO_READER[filetype]

    # Open and return a filehandle suitable for reading this archive.
    return reader(filename)


@type_check
def write_bytes(filename: str, is_overwritable: bool = False) -> BufferedIOBase:
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
    is_overwritable : optional[bool]
        `True` if overwriting this file when this file already exists _or_
        `False` if raising an exception when this file already exists. Defaults
        to `False` for safety.

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
    from betse.util.path import dirs, paths, pathnames

    # Raise an exception unless this filename has a supported archive filetype.
    die_unless_filetype(filename)

    # If this file is *NOT* overwritable, raise an exception if this path
    # already exists.
    if not is_overwritable:
        paths.die_if_path(filename)

    # Create the parent directory of this file if needed.
    dirs.make_parent_unless_dir(filename)

    # Filetype of this pathname. By the above validation, this filetype is
    # guaranteed to both exist and be a supported archive filetype. Hence, no
    # additional validation is required.
    filetype = pathnames.get_filetype_undotted(filename)

    # Callable writing archives of this filetype.
    writer = _ARCHIVE_FILETYPE_TO_WRITER[filetype]

    # Open and return a filehandle suitable for writing this archive.
    return writer(filename, is_overwritable=is_overwritable)

# ....................{ IO ~ bz2                           }....................
def _read_bytes_bz2(filename: str) -> BufferedIOBase:
    '''
    Open and return a filehandle suitable for reading the binary bzip-archived
    file with the passed filename.

    This private function is intended to be called *only* by the public
    :func:`reading_bytes` function.
    '''

    # This optional stdlib module is guaranteed to exist and hence be safely
    # importable here, due to the above die_unless_filetype() call.
    from bz2 import BZ2File

    # Open and return a filehandle suitable for reading this file.
    return BZ2File(filename, mode='rb')


def _write_bytes_bz2(filename: str, is_overwritable: bool) -> BufferedIOBase:
    '''
    Open and return a filehandle suitable for writing the binary bzip-archived
    file with the passed filename.

    This private function is intended to be called *only* by the public
    :func:`writing_bytes` function.
    '''

    # Avoid circular import dependencies.
    from betse.util.io import iofiles

    # This optional stdlib module is guaranteed to exist and hence be safely
    # importable here, due to the above die_unless_filetype() call.
    from bz2 import BZ2File

    # Open and return a filehandle suitable for e(x)clusively writing this file.
    return BZ2File(filename, mode=iofiles.get_mode_write_bytes(is_overwritable))

# ....................{ IO ~ gz                            }....................
def _read_bytes_gz(filename: str) -> BufferedIOBase:
    '''
    Open and return a filehandle suitable for reading the binary gzip-archived
    file with the passed filename.

    This private function is intended to be called *only* by the public
    :func:`reading_bytes` function.
    '''

    # This optional stdlib module is guaranteed to exist and hence be safely
    # importable here, due to the above die_unless_filetype() call.
    from gzip import GzipFile

    # Open and return a filehandle suitable for reading this file.
    return GzipFile(filename, mode='rb')


def _write_bytes_gz(filename: str, is_overwritable: bool) -> BufferedIOBase:
    '''
    Open and return a filehandle suitable for writing the binary gzip-archived
    file with the passed filename.

    This private function is intended to be called *only* by the public
    :func:`writing_bytes` function.
    '''

    # Avoid circular import dependencies.
    from betse.util.io import iofiles

    # This optional stdlib module is guaranteed to exist and hence be safely
    # importable here, due to the above die_unless_filetype() call.
    from gzip import GzipFile

    # Open and return a filehandle suitable for e(x)clusively writing this file.
    return GzipFile(filename, mode=iofiles.get_mode_write_bytes(is_overwritable))

# ....................{ IO ~ xz                            }....................
def _read_bytes_xz(filename: str) -> BufferedIOBase:
    '''
    Open and return a filehandle suitable for reading the binary LZMA-archived
    file with the passed filename.

    This private function is intended to be called *only* by the public
    :func:`reading_bytes` function.
    '''

    # This optional stdlib module is guaranteed to exist and hence be safely
    # importable here, due to the above die_unless_filetype() call.
    from lzma import LZMAFile

    # Open and return a filehandle suitable for reading this file.
    return LZMAFile(filename, mode='rb')


def _write_bytes_xz(filename: str, is_overwritable: bool) -> BufferedIOBase:
    '''
    Open and return a filehandle suitable for writing the binary LZMA-archived
    file with the passed filename.

    This private function is intended to be called *only* by the public
    :func:`writing_bytes` function.
    '''

    # Avoid circular import dependencies.
    from betse.util.io import iofiles

    # This optional stdlib module is guaranteed to exist and hence be safely
    # importable here, due to the above die_unless_filetype() call.
    from lzma import LZMAFile

    # Open and return a filehandle suitable for e(x)clusively writing this file.
    return LZMAFile(filename, mode=iofiles.get_mode_write_bytes(is_overwritable))

# ....................{ CONSTANTS ~ private                }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: When adding a new key-value pair or removing an existing key-value
# pair from the dictionaries defined below, all testing parametrizations
# exercising these dictionaries (e.g., in the
# "betse_test.unit.path.test_archive" submodule) *MUST* be updated accordingly.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

_ARCHIVE_FILETYPE_TO_MODULE_NAME = {
    'bz2': 'bz2',
    'gz':  'gzip',
    'xz':  'lzma',
}
'''
Dictionary mapping from each archive filetype supported by this submodule to
the fully-qualified name of the optional stdlib module required to write and
read archives of this filetype.
'''


_ARCHIVE_FILETYPE_TO_READER = {
    'bz2': _read_bytes_bz2,
    'gz':  _read_bytes_gz,
    'xz':  _read_bytes_xz,
}
'''
Dictionary mapping from each archive filetype supported by this submodule to
the private function of this submodule reading archives of this filetype.
'''


_ARCHIVE_FILETYPE_TO_WRITER = {
    'bz2': _write_bytes_bz2,
    'gz':  _write_bytes_gz,
    'xz':  _write_bytes_xz,
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
    from betse.util.type import modules

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
