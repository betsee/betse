#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level non-directory file facilities.

This module is named `files` rather than `file` to avoid conflict with the stock
`files` class.
'''

# ....................{ IMPORTS                            }....................
import os, re, shutil
from betse.exceptions import BetseFileException
from betse.util.io.log import logs
from betse.util.type.types import type_check, SequenceTypes
from io import BufferedIOBase, TextIOWrapper
from os import path

# ....................{ EXCEPTIONS ~ unless                }....................
def die_unless_file(pathname: str) -> None:
    '''
    Raise an exception unless the passed path is a non-directory file *after*
    following symbolic links.

    See Also
    ----------
    :func:`is_file`
        Further details.
    '''

    if not is_file(pathname):
        raise BetseFileException(
            'File "{}" not found or unreadable.'.format(pathname))

# ....................{ EXCEPTIONS ~ if                    }....................
def die_if_file(pathname: str) -> None:
    '''
    Raise an exception if the passed path is a non-directory file *after*
    following symbolic links.

    See Also
    ----------
    :func:`is_file`
        Further details.
    '''

    if is_file(pathname):
        raise BetseFileException('File "{}" already exists.'.format(pathname))


def die_if_special(pathname: str) -> None:
    '''
    Raise an exception if the passed path is a special file.

    See Also
    ----------
    :func:`is_special`
        Further details.
    '''

    if is_special(pathname):
        # Avoid circular import dependencies.
        from betse.util.path import paths

        raise BetseFileException(
            'File "{}" already an existing {}.'.format(
                pathname, paths.get_type_label(pathname)))

# ....................{ TESTERS                            }....................
def is_file(pathname: str) -> bool:
    '''
    `True` only if the passed path is a non-directory file _after_ following
    symbolic links.

    This function does _not_ raise an exception if this path does not exist.

    Versus `path.isfile()`
    ----------
    This function intrinsically differs from the standard :func:`path.isfile`
    function. While the latter returns `True` only for non-special files and
    hence `False` for all non-directory special files (e.g., device nodes,
    sockets), this function returns `True` for _all_ non-directory files
    regardless of whether these files are special or not.

    **Why?** Because this function complies with POSIX semantics, whereas
    :func:`path.isfile` does _not_. The specialness of non-directory files is
    usually irrelevant; in general, it only matters whether these files are
    directories or not. For example, the external command `rm` removes only
    non-directory files (regardless of specialness) while the external command
    `rmdir` removes only empty directories.
    '''

    # Avoid circular import dependencies.
    from betse.util.path import dirs, paths

    # This path is a file if this path both exists and is *NOT* a directory.
    return paths.is_path(pathname) and not dirs.is_dir(pathname)


def is_file_executable(pathname: str) -> bool:
    '''
    `True` only if the passed path is an **executable non-directory file**
    (i.e., file with the execute bit enabled) _after_ following symbolic links.

    This function does _not_ raise an exception if this path does not exist.
    '''

    # This path is an executable file if this path is an existing file with the
    # executable bit enabled.
    return is_file(pathname) and os.access(pathname, os.X_OK)


#FIXME: Given that this function resides in the "files" submodule, shouldn't
#this function return False rather than True when passed a directory pathname?
def is_special(pathname: str) -> bool:
    '''
    `True` only if the passed path is an existing **special file** (e.g.,
    directory, device node, socket, symbolic link).

    This function does _not_ raise an exception if this path does not exist.
    '''

    # Avoid circular import dependencies.
    from betse.util.path import paths

    # True if this path exists and...
    return paths.is_path(pathname) and (
        # ...is either a symbolic link *OR* neither a regular file nor symbolic
        # link to such a file. In the latter case, predicate logic guarantees
        # this file to *NOT* be a symbolic link, thus reducing this test to:
        # "...is either a symbolic link *OR* not a regular file."
        is_symlink(pathname) or not path.isfile(pathname))

# ....................{ TESTERS ~ symlink                  }....................
@type_check
def is_symlink(pathname: str) -> bool:
    '''
    `True` only if the passed path is an existing symbolic link _before_
    following symbolic links.

    This function does _not_ raise an exception if this path does not exist.
    '''

    return path.islink(pathname)


@type_check
def is_symlink_valid(pathname: str) -> bool:
    '''
    `True` only if the passed path is an existing **non-dangling symbolic link**
    (i.e., symbolic link whose target also exists) _before_ following symbolic
    links.

    This function does _not_ raise an exception if this path does not exist.
    '''

    # Call path.exists() rather than path.lexists(), as the latter returns True
    # for dangling symbolic links.
    #
    # This is why human-readable function names is a good thing, people.
    return is_symlink(pathname) and path.exists(pathname)

# ....................{ COPIERS                            }....................
@type_check
def get_chars(filename: str, encoding: str = 'utf-8') -> str:
    '''
    String of all characters contained in the plaintext file with the passed
    filename encoded with the passed encoding.

    Parameters
    ----------
    filename : str
        Relative or absolute path of the plaintext text to be read.
    encoding : optional[str]
        Name of the encoding to be used. Defaults to UTF-8.

    Returns
    ----------
    str
        String of all characters decoded from this file's byte content.
    '''

    with read_chars(filename=filename, encoding=encoding) as text_file:
        return text_file.read()


@type_check
def get_mode_write_bytes(is_overwritable: bool = False) -> str:
    '''
    Mode string suitable for opening a file handle for byte-oriented writing via
    the `mode` parameter to the :func:`open` builtin.

    This low-level utility function is intended to be called _only_ by
    higher-level utility functions (e.g., :func:`write_bytes`).

    Parameters
    ----------
    is_overwritable : optional[bool]
        `True` if overwriting this file when this file already exists _or_
        `False` if raising an exception when this file already exists. Defaults
        to `False` for safety.

    Returns
    ----------
    str
        If `is_overwritable` is:
        * `True`, this is `"xb"`, raising exceptions when attempting to write
          files that already exist with this mode.
        * `False`, this is `"wb"`, silently overwriting such files.
    '''

    return 'wb' if is_overwritable else 'xb'

# ....................{ COPIERS                            }....................
@type_check
def copy(filename_source: str, filename_target: str) -> None:
    '''
    Copy the passed source file to the passed target file or directory.

    If the source file is a symbolic link, this link (rather than its
    transitive target) will be copied and hence preserved.

    The target file will be copied in a manner maximally preserving metadata
    (e.g., owner, group, permissions, times, extended file system attributes).
    If the target file is a directory, the basename of the source file will be
    appended to this directory -- much like the standard `cp` POSIX command.

    If either the source file does not exist _or_ the target file already
    exists, an exception will be raised.
    '''

    # Avoid circular import dependencies.
    from betse.util.path import dirs, paths

    # Log this copy.
    logs.log_debug(
        'Copying file "%s" to "%s".', filename_source, filename_target)

    # Raise an exception unless the source file exists.
    die_unless_file(filename_source)

    # If the target file is a directory, append the basename of the passed
    # source file to this directory -- much like the "cp" POSIX command.
    if dirs.is_dir(filename_target):
        filename_target = paths.join(
            filename_target, paths.get_basename(filename_source))

    # Raise an exception if the target file already exists.
    paths.die_if_path(filename_target)

    # Perform this copy in a manner preserving metadata and symbolic links.
    shutil.copy2(filename_source, filename_target, follow_symlinks=False)

# ....................{ REMOVERS                           }....................
@type_check
def remove_if_found(filename: str) -> None:
    '''
    Remove the passed non-directory file if this file currently exists.

    If this file does _not_ currently exist, this function reduces to a noop.
    For safety, this function removes this file atomically; in particular, this
    file's existence is _not_ explicitly tested for.
    '''

    # Log this removal if the subsequent removal attempt is likely to actually
    # remove a file. Due to race conditions with other threads and processes,
    # this file could be removed after this test succeeds but before the removal
    # is performed. Since this is largely ignorable, the worst case is an
    # extraneous log message.
    if is_file(filename):
        logs.log_debug('Removing file: %s', filename)

    # Remove this file atomically. To avoid race conditions with other threads
    # and processes, this operation is *NOT* embedded in an explicit test for
    # file existence. Instead, the Pythonic Way is embraced.
    try:
        os.remove(filename)
    # If this file does *NOT* exist, ignore this exception.
    except FileNotFoundError:
        pass


@type_check
def remove(filename: str) -> None:
    '''
    Remove the passed non-directory file.
    '''

    # Log this removal.
    logs.log_debug('Removing file: %s', filename)

    # Raise an exception unless this file exists.
    die_unless_file(filename)

    # Remove this file. Note that the os.remove() and os.unlink() functions are
    # identical. (That was a tad silly, Guido.)
    os.remove(filename)

# ....................{ READERS                            }....................
@type_check
def read_bytes(filename: str) -> BufferedIOBase:
    '''
    Open and return a filehandle suitable for reading the binary file with the
    passed filename, transparently decompressing this file if the filetype of
    this filename is that of a supported archive format.

    This function returns a :class:`file`-like object suitable for use wherever
    the :func:`open` builtin is callable (e.g., in `with` statements).

    Parameters
    ----------
    filename : str
        Relative or absolute path of the binary file to be read. If this
        filename is suffixed by a supported archive filetype (i.e., if the
        :func:`betse.util.path.archives.is_filetype` function returns `True`
        for this filename), the returned filehandle automatically reads the
        decompressed rather than compressed byte contents of this file.

    Returns
    ----------
    BufferedIOBase
        `file`-like object encapsulating this opened file.

    Raises
    ----------
    BetseFileException
        If this filename is _not_ that of an existing file.
    '''

    # Avoid circular import dependencies.
    from betse.util.path import archives

    # Log this I/O operation.
    logs.log_debug('Reading bytes: %s', filename)

    # Raise an exception unless this file exists.
    die_unless_file(filename)

    # If this file is compressed, open and return a file handle reading
    # decompressed bytes from this file.
    if archives.is_filetype(filename):
        return archives.read_bytes(filename)
    # Else, this file is uncompressed. Open and return a typical file handle
    # reading bytes from this file.
    else:
        return open(filename, mode='rb')


@type_check
def read_chars(filename: str, encoding: str = 'utf-8') -> TextIOWrapper:
    '''
    Open and return a filehandle suitable for reading the plaintext file with
    the passed filename encoded with the passed encoding.

    This function returns a :class:`file`-like object suitable for use wherever
    the :func:`open` builtin is callable (e.g., in `with` statements).

    Parameters
    ----------
    filename : str
        Relative or absolute path of the plaintext text to be read.
    encoding : optional[str]
        Name of the encoding to be used. Defaults to UTF-8.

    Returns
    ----------
    TextIOWrapper
        `file`-like object encapsulating the opened file.
    '''

    # Log this I/O operation.
    logs.log_debug('Reading chars: %s', filename)

    # Raise an exception unless this file exists.
    die_unless_file(filename)

    # Open this file.
    return open(filename, mode='rt', encoding=encoding)

# ....................{ WRITERS                            }....................
@type_check
def write_bytes(filename: str, is_overwritable: bool = False) -> BufferedIOBase:
    '''
    Open and return a filehandle suitable for writing the binary file with the
    passed filename, transparently compressing this file if the filetype of
    this filename is that of a supported archive format.

    This function returns a :class:`file`-like object suitable for use wherever
    the :func:`open` builtin is callable (e.g., in `with` statements).

    Parameters
    ----------
    filename : str
        Relative or absolute path of the binary file to be written. If this
        filename is suffixed by a supported archive filetype (i.e., if the
        :func:`betse.util.path.archives.is_filetype` function returns `True`
        for this filename), the returned filehandle automatically writes the
        compressed rather than uncompressed byte contents of this file.
    is_overwritable : optional[bool]
        `True` if overwriting this file when this file already exists _or_
        `False` if raising an exception when this file already exists. Defaults
        to `False` for safety.

    Returns
    ----------
    BufferedIOBase
        `file`-like object encapsulating this opened file.
    '''

    # Avoid circular import dependencies.
    from betse.util.path import archives, dirs, paths

    # Log this I/O operation.
    logs.log_debug('Writing bytes: %s', filename)

    # If this file is *NOT* overwritable, raise an exception if this path
    # already exists.
    if not is_overwritable:
        paths.die_if_path(filename)

    # Create the parent directory of this file if needed.
    dirs.make_parent_unless_dir(filename)

    # If this file is compressed, open and return a file handle writing
    # compressed bytes to this file.
    if archives.is_filetype(filename):
        return archives.write_bytes(filename, is_overwritable=is_overwritable)
    # Else, this file is uncompressed.
    else:
        # Mode with which to open this file for byte-oriented writing.
        mode = get_mode_write_bytes(is_overwritable)

        # Open and return a file handle writing uncompressed bytes to this file.
        return open(filename, mode=mode)


#FIXME: Add the "is_overwritable" parameter, implemented similarly to the same
#parameter accepted by the write_bytes() function.
@type_check
def write_chars(filename: str, encoding: str = 'utf-8') -> TextIOWrapper:
    '''
    Open and return a filehandle suitable for writing the plaintext file with
    the passed filename encoded with the passed encoding.

    This function returns a :class:`file`-like object suitable for use wherever
    the :func:`open` builtin is callable (e.g., in `with` statements).

    Parameters
    ----------
    filename : str
        Relative or absolute path of the plaintext file to be written.
    encoding : optional[str]
        Name of the encoding to be used. Defaults to UTF-8.

    Returns
    ----------
    TextIOWrapper
        `file`-like object encapsulating this opened file.
    '''

    # Avoid circular import dependencies.
    from betse.util.path import dirs, paths

    # Log this I/O operation.
    logs.log_debug('Writing chars: %s', filename)

    # Raise an exception if this path already exists.
    paths.die_if_path(filename)

    # Create the parent directory of this file if needed.
    dirs.make_parent_unless_dir(filename)

    # Open this file.
    return open(filename, mode='xt', encoding=encoding)

# ....................{ REPLACERS                          }....................
def replace_substrs_inplace(
    filename: str, replacements: SequenceTypes, **kwargs) -> None:
    '''
    Replace all substrings in the passed non-directory file matching the passed
    regular expressions with the corresponding passed substitutions.

    See Also
    ----------
    :func:`replace_substrs`
        For further details.
    '''

    replace_substrs(filename, filename, replacements, **kwargs)


@type_check
def replace_substrs(
    filename_source: str,
    filename_target: str,
    replacements: SequenceTypes,
    **kwargs
) -> None:
    '''
    Write the passed target non-directory file with the result of replacing all
    substrings in the passed source non-directory file matching the passed
    regular expressions with the corresponding passed substitutions.

    This function implements the equivalent of the `sed` line processor in a
    pure-Python manner requiring no external commands or additional
    dependencies. This is a good thing.

    Such source file will _not_ be changed. Such target file will be written in
    an atomic manner maximally preserving source metadata (e.g., owner, group,
    permissions, times, extended file system attributes). If either the source
    file does not exist _or_ the target file already exists, an exception will
    be raised.

    Such regular expressions may be either strings *or* instances of the
    `Pattern` class (i.e., compiled regular expression object), while such
    substitutions may be either strings _or_ functions. See the "Arguments"
    section below for how this function expects such objects to be passed.

    This function accepts the same optional keyword arguments as `re.sub()`.

    Arguments
    ----------
    filename_source : str
        Absolute path of the source filename to be read.
    filename_target : str
        Absolute path of the target filename to be written.
    replacements : Sequence
        Non-string sequence (e.g., list, tuple) of non-string sequences of
        length 2 (i.e., pairs), whose:
        * First element is a regular expression matching all substrings in the
          source file to be replaced.
        * Second element is the substring in the target file to replace all
          substrings in the source file matching that regular expression.
    '''

    # Avoid circular import dependencies.
    from betse.util.path import temps

    # Log this substitution.
    if filename_source == filename_target:
        logs.log_debug(
            'Munging file "%s" in-place.', filename_source)
    else:
        logs.log_debug(
            'Munging file "%s" to "%s".', filename_source, filename_target)

    # Raise an exception unless the source file exists.
    die_unless_file(filename_source)

    # Raise an exception if the target file already exists.
    die_if_file(filename_target)

    # For efficiency, replace all passed uncompiled with compiled regular
    # expressions via a tuple comprehension.
    replacements = tuple(
        (re.compile(regex), substitution)
        for (regex, substitution) in replacements
    )

    #FIXME: This functionality is probably quite useful, where at least one
    #matching line is absolutely expected. Consider formalizing into a passed
    #argument, as lackluster time permits.
    # is_line_matches = False

    # For atomicity, incrementally write to a temporary file rather than the
    # desired target file *BEFORE* moving the former to the latter. This
    # obscures such writes from other threads and/or processes, avoiding
    # potential race conditions elsewhere.
    with temps.write_chars() as file_target_temp:
        with read_chars(filename_source) as file_source:
            # For each line of the source file...
            for line in file_source:
                # For each passed regular expression and corresponding
                # substitution, replace all substrings in this line matching
                # that regular expression with that substitution.
                for (regex, substitution) in replacements:
                    # if regex.search(line, **kwargs) is not None:
                    #     is_line_matches = True
                    #     loggers.log_info('Line "%s" matches!', line)
                    line = regex.sub(substitution, line, **kwargs)

                # Append such line to the temporary file.
                file_target_temp.write(line)

            # if not is_line_matches:
            #     raise BetseFileException('No line matches!')

    # Copy all metadata (e.g., permissions) from the source to temporary
    # file *BEFORE* moving the latter, avoiding potential race conditions
    # and security vulnerabilities elsewhere.
    shutil.copystat(filename_source, file_target_temp.name)

    # Move the temporary to the target file.
    shutil.move(file_target_temp.name, filename_target)
