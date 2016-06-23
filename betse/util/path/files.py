#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level non-directory file facilities.

This module is named `files` rather than `file` to avoid conflict with the stock
`files` class.
'''

# ....................{ IMPORTS                            }....................
import os, re, shutil, tempfile
from betse.exceptions import BetseExceptionFile
from betse.util.io.log import logs
from betse.util.type import types
from betse.util.type.types import type_check
from collections.abc import Sequence
from os import path

# ....................{ EXCEPTIONS ~ unless                }....................
def die_unless_file(pathname: str) -> None:
    '''
    Raise an exception unless the passed non-directory file exists _after_
    following symbolic links.

    See Also
    ----------
    `is_file()`
        For further details.
    '''
    if not is_file(pathname):
        raise BetseExceptionFile(
            'File "{}" not found or unreadable.'.format(pathname))

# ....................{ EXCEPTIONS ~ if                    }....................
def die_if_file(pathname: str) -> None:
    '''
    Raise an exception if the passed non-directory file exists _after_ following
    symbolic links.

    See Also
    ----------
    `is_file()`
        For further details.
    '''
    if is_file(pathname):
        raise BetseExceptionFile('File "{}" already exists.'.format(pathname))


def die_if_special(pathname: str) -> None:
    '''
    Raise an exception if the passed path is an existing special path.

    See Also
    ----------
    `is_special()`
        For further details.
    '''
    if is_special(pathname):
        # Avoid circular import dependencies.
        from betse.util.path import paths
        raise BetseExceptionFile(
            'File "{}" already an existing {}.'.format(
                pathname, paths.get_type_label(pathname)))

# ....................{ TESTERS                            }....................
def is_file(pathname: str) -> bool:
    '''
    `True` if the passed path is an existing non-directory file exists _after_
    following symbolic links.

    Versus `path.isfile()`
    ----------
    This function intrinsically differs from the standard `path.isfile()`
    function. While the latter returns `True` only for non-special files and
    hence `False` for all non-directory special files (e.g., device nodes,
    sockets), this function returns `True` for _all_ non-directory files
    regardless of whether such files are special or not.

    **Why?** Because this function complies with POSIX semantics, whereas
    `path.isfile()` does _not_. The specialness of non-directory files is
    usually irrelevant; in general, it only matters whether such files are
    directories or not. For example, the external command `rm` removes only
    non-directory files (regardless of specialness) while the external command
    `rmdir` removes only empty directories.
    '''
    # Avoid circular import dependencies.
    from betse.util.path import dirs, paths

    # True if such path exists and is *NOT* a directory.
    return paths.is_path(pathname) and not dirs.is_dir(pathname)


def is_special(pathname: str) -> bool:
    '''
    `True` if the passed path is an existing **special file** (e.g., directory,
    device node, socket, symbolic link).
    '''
    # Avoid circular import dependencies.
    from betse.util.path import paths

    # True if such path exists and...
    return paths.is_path(pathname) and (
        # ...is either a symbolic link *OR* neither a regular file nor symbolic
        # link to such a file. In the latter case, predicate logic guarantees
        # such file to *NOT* be a symbolic link, thus reducing this test to:
        # "...is either a symbolic link *OR* not a regular file."
        is_symlink(pathname) or not path.isfile(pathname))


def is_symlink(pathname: str) -> bool:
    '''
    `True` if the passed path is an existing symbolic link.
    '''
    assert types.is_str_nonempty(pathname),\
        types.assert_not_str_nonempty(pathname, 'pathname')
    return path.islink(pathname)

# ....................{ COPIERS                            }....................
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
    assert types.is_str_nonempty(filename_source), (
        types.assert_not_str_nonempty(filename_source, 'source filename'))
    assert types.is_str_nonempty(filename_target), (
        types.assert_not_str_nonempty(filename_target, 'target filename'))

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
def remove(filename: str) -> None:
    '''
    Remove the passed non-directory file.
    '''
    assert types.is_str_nonempty(filename), (
        types.assert_not_str_nonempty(filename, 'filename'))

    # Log this removal.
    logs.log_debug('Removing file "%s".', filename)

    # Raise an exception unless such this exists.
    die_unless_file(filename)

    # Remove this file. Note that the os.remove() and os.unlink() functions are
    # identical. (That was silly, Guido.)
    os.remove(filename)


def remove_if_found(filename: str) -> None:
    '''
    Remove the passed non-directory file if this file currently exists.

    If this file does _not_ currently exist, this function reduces to a noop.
    For safety, this function removes this file atomically; in particular, this
    file's existence is _not_ explicitly tested for.
    '''
    assert types.is_str_nonempty(filename), (
        types.assert_not_str_nonempty(filename, 'filename'))

    # Log this removal if the subsequent removal attempt is likely to actually
    # remove a file. Due to race conditions with other processes, this file
    # could be removed after this test succeeds but before the removal is
    # performed. Since this is largely ignorable, the worst case is an
    # extraneous log message.
    if is_file(filename):
        logs.log_debug('Removing file "%s".', filename)

    # Remove this file atomically. To avoid race conditions with other
    # processes, do *NOT* embed this operation in an explicit test for file
    # existence. Instead, adopt the Pythonic Way.
    try:
        os.remove(filename)
    # If this file does *NOT* exist, ignore this exception.
    except FileNotFoundError:
        pass

# ....................{ OPENERS                            }....................
def open_for_text_reading(filename: str) -> 'file':
    '''
    Open and return the passed file for line-oriented reading.

    This function returns a `file`-like object, suitable for use wherever the
    builtin `open()` would otherwise be called (e.g., in `with` statements).

    Returns
    ----------
    file
        `file`-like object encapsulating the opened file.
    '''
    assert types.is_str_nonempty(filename), (
        types.assert_not_str_nonempty(filename, 'filename'))

    # Raise an exception unless this file exists.
    die_unless_file(filename)

    # Open this file.
    return open(filename, mode='rt')


def open_for_text_writing(filename: str, encoding='utf-8'):
    '''
    Open and return the passed file for line-oriented writing.

    This function returns a `file`-like object, suitable for use wherever the
    builtin `open()` would otherwise be called (e.g., in `with` statements).

    Parameters
    ----------
    encoding : optional[str]
        Optional encoding to be used. Defaults to UTF-8.

    Returns
    ----------
    file
        `file`-like object encapsulating the opened file.
    '''
    assert types.is_str_nonempty(filename), (
        types.assert_not_str_nonempty(filename, 'filename'))

    # Avoid circular import dependencies.
    from betse.util.path import dirs

    # Raise an exception if this path is an existing special file. Such files
    # must *NEVER* be opened for line-oriented writing.
    die_if_special(filename)

    # Create the parent directory of this file if needed.
    dirs.make_parent_unless_dir(filename)

    # Open this file.
    return open(filename, mode='wt', encoding=encoding)

# ....................{ OPENERS ~ temporary                }....................
def open_for_text_writing_temporary(encoding='utf-8'):
    '''
    Open and return a temporary named file for line-oriented writing.

    This function returns a `file`-like object, suitable for use wherever the
    builtin `open()` would otherwise be called (e.g., in `with` statements).
    The filename for the underlying file is accessible in the typical way for
    named `file`-like objects (namely, via the `name` attribute of such object).
    The underlying file will _not_ be implicitly deleted when such object is
    closed, which typically defeats the purpose of creating a temporary named
    file in the first place.

    Text Lines
    ----------
    Temporary files are _not_ reliably openable for line-oriented writing in a
    cross-platform manner. Hence, callers are discouraged from calling this
    function except where absolutely required (e.g., some or all underlying
    encodings are unknown). Instead, text lines may be manually written in a
    byte-oriented manner by encoding such lines under a known encoding _and_
    appending such lines by the newline delimiter specific to the current
    platform: e.g.,

        >>> from betse.util.path imports files
        >>> import os
        >>> tempfile = files.open_for_byte_writing_temporary(encoding='utf-8')
        >>> tempfile.write(bytes('Eaarth' + os.linesep))

    Parameters
    ----------
    encoding : optional[str]
        Optional encoding to be used. Defaults to UTF-8.

    Returns
    ----------
    file
        `file`-like object encapsulating the opened file.
    '''

    return tempfile.NamedTemporaryFile(
        mode='w+', delete=False, encoding=encoding)


def open_for_byte_writing_temporary(encoding='utf-8'):
    '''
    Open and return a temporary named file for byte-oriented writing.

    Parameters
    ----------
    encoding : optional[str]
        Optional encoding to be used. Defaults to UTF-8.

    See Also
    ----------
    `open_for_text_writing_temporary()`
        For further details.

    Returns
    ----------
    file
        `file`-like object encapsulating the opened file.
    '''

    return tempfile.NamedTemporaryFile(delete=False, encoding=encoding)

# ....................{ OPENERS ~ temporary                }....................
def substitute_substrings_inplace(
    filename: str, substitutions, **kwargs) -> None:
    '''
    Replace all substrings in the passed non-directory file matching the passed
    regular expressions with the corresponding passed substitutions.

    See Also
    ----------
    `substitute_substrings()`
        For further details.
    '''

    substitute_substrings(
        filename, filename, substitutions, **kwargs)


@type_check
def substitute_substrings(
    filename_source: str,
    filename_target: str,
    substitutions: Sequence,
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
    substitutions : Sequence
        Non-string sequence (e.g., list, tuple) of non-string sequences of
        length 2 (i.e., pairs), whose first element is a regular expression and
        whose second element is the substitution to be performed for all
        substrings in the source file matching that regular expression.
    '''

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
    substitutions = tuple(
        (re.compile(regex), substitution)
        for (regex, substitution) in substitutions
    )

    #FIXME: This functionality is probably quite useful, where at least one
    #matching line is absolutely expected. Consider formalizing into a passed
    #argument, as lackluster time permits.
    # is_line_matches = False

    # For atomicity, incrementally write to a temporary file rather than the
    # desired target file *BEFORE* moving the former to the latter. This
    # obscures such writes from other threads and/or processes, avoiding
    # potential race conditions elsewhere.
    with open_for_text_writing_temporary() as file_target_temp:
        with open_for_text_reading(filename_source) as file_source:
            # For each line of the source file...
            for line in file_source:
                # For each passed regular expression and corresponding
                # substitution, replace all substrings in this line matching
                # that regular expression with that substitution.
                for (regex, substitution) in substitutions:
                    # if regex.search(line, **kwargs) is not None:
                    #     is_line_matches = True
                    #     loggers.log_info('Line "%s" matches!', line)
                    line = regex.sub(substitution, line, **kwargs)

                # Append such line to the temporary file.
                file_target_temp.write(line)

            # if not is_line_matches:
            #     raise BetseExceptionFile('No line matches!')

    # Copy all metadata (e.g., permissions) from the source to temporary
    # file *BEFORE* moving the latter, avoiding potential race conditions
    # and security vulnerabilities elsewhere.
    shutil.copystat(filename_source, file_target_temp.name)

    # Move the temporary to the target file.
    shutil.move(file_target_temp.name, filename_target)
