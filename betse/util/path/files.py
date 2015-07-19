#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level non-directory file facilities.

This module is named `files` rather than `file` to avoid conflict with the stock
`files` class.
'''

# ....................{ IMPORTS                            }....................
from betse.exceptions import BetseExceptionFile
from betse.util.io import loggers
from os import path
import os, re, shutil, tempfile

# ....................{ EXCEPTIONS ~ unless                }....................
def die_unless_file(pathname: str) -> None:
    '''
    Raise an exception unless the passed non-directory file exists *after*
    following symbolic links.
    '''
    if not is_file(pathname):
        raise BetseExceptionFile(
            'File "{}" not found or unreadable.'.format(pathname))

# ....................{ EXCEPTIONS ~ if                    }....................
def die_if_file(pathname: str) -> None:
    '''
    Raise an exception if the passed non-directory file exists *after* following
    symbolic links.
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
    `True` if the passed path is an existing non-directory file exists *after*
    following symbolic links.

    Versus `path.isfile()`
    ----------
    This function intrinsically differs from the standard `path.isfile()`
    function. While the latter returns `True` only for non-special files and
    hence `False` for all non-directory special files (e.g., device nodes,
    sockets), this function returns `True` for *all* non-directory files
    regardless of whether such files are special or not.

    Why? Because this function complies with POSIX semantics, while
    `path.isfile()` does *not*. The specialness of non-directory files is
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
    assert isinstance(pathname, str), '"{}" not a string.'.format(pathname)
    assert len(pathname), 'Pathname empty.'
    return path.islink(pathname)

# ....................{ COPIERS                            }....................
def copy(filename_source: str, filename_target: str) -> None:
    '''
    Copy the passed source to target non-directory file.

    Such file will be copied in a manner maximally preserving metadata (e.g.,
    owner, group, permissions, times, extended file system attributes).
    Likewise, if such source file is a symbolic link, such link (rather than its
    transitive target) will be copied and hence preserved.

    If either the source file does not exist *or* the target file already
    exists, an exception will be raised.
    '''
    assert isinstance(filename_source, str),\
        '"{}" not a string.'.format(filename_source)
    assert isinstance(filename_target, str),\
        '"{}" not a string.'.format(filename_target)
    assert len(filename_source), 'Source filename empty.'
    assert len(filename_target), 'Target filename empty.'

    # Log such copy.
    loggers.log_info(
        'Copying file "%s" to "%s".', filename_source, filename_target)

    # Raise an exception unless the source file exists.
    die_unless_file(filename_source)

    # Raise an exception if the target file already exists.
    die_if_file(filename_target)

    # Perform such copy in a manner preserving metadata and symbolic links.
    shutil.copy2(filename_source, filename_target, follow_symlinks = False)

# ....................{ REMOVERS                           }....................
def remove(filename: str) -> None:
    '''
    Remove the passed non-directory file.
    '''
    assert isinstance(filename, str), '"{}" not a string.'.format(filename)
    assert len(filename), 'Filename empty.'

    # Log such removal.
    loggers.log_info('Removing file "%s".', filename)

    # Raise an exception unless such file exists.
    die_unless_file(filename)

    # Remove such file. Note that the os.remove() and os.unlink() functions are
    # identical. (That was silly, Guido.)
    os.remove(filename)

# ....................{ OPENERS                            }....................
def open_for_text_reading(filename: str):
    '''
    Open the passed file for line-oriented reading.

    This function returns a `file`-like object, suitable for use wherever the
    builtin `open()` would otherwise be called (e.g., in `with` statements).
    '''
    assert isinstance(filename, str), '"{}" not a string.'.format(filename)
    assert len(filename), 'Filename empty.'

    # Raise an exception unless such file exists.
    die_unless_file(filename)

    # Open such file.
    return open(filename, mode = 'rt')

def open_for_text_writing(filename: str):
    '''
    Open the passed file for line-oriented writing.

    This function returns a `file`-like object, suitable for use wherever the
    builtin `open()` would otherwise be called (e.g., in `with` statements).
    '''
    assert isinstance(filename, str), '"{}" not a string.'.format(filename)
    assert len(filename), 'Filename empty.'

    # Avoid circular import dependencies.
    from betse.util.path import dirs

    # Raise an exception if such path is an existing special file. Such files
    # must *NEVER* be opened for line-oriented writing.
    die_if_special(filename)

    # Create the parent directory of such file if needed.
    dirs.make_parent_unless_dir(filename)

    # Open such file.
    return open(filename, mode = 'wt')

# ....................{ OPENERS ~ temporary                }....................
def open_for_text_writing_temporary():
    '''
    Open a temporary named file for line-oriented writing.

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
        >>> tempfile = files.open_for_byte_writing_temporary()
        >>> tempfile.write(bytes('Eaarth' + os.linesep, encoding = 'utf-8'))
    '''
    return tempfile.NamedTemporaryFile(mode='w+', delete=False)

def open_for_byte_writing_temporary():
    '''
    Open a temporary named file for byte-oriented writing.

    See Also
    ----------
    `open_for_text_writing_temporary()`
        For further details.
    '''
    return tempfile.NamedTemporaryFile(delete=False)

# ....................{ OPENERS ~ temporary                }....................
def substitute_strings(
    filename_source: str,
    filename_target: str,
    regex,
    substitution,
    **kwargs
) -> None:
    '''
    Write the passed target non-directory file with the result of substituting
    all substrings in the passed source non-directory file matching the passed
    regular expression with the passed substitution.

    This function implements the equivalent of the `sed` line processor in a
    pure-Python manner requiring no external commands or additional
    dependencies. This is a good thing.

    Such source file will _not_ be changed. Such target file will be written in
    an atomic manner maximally preserving source metadata (e.g., owner, group,
    permissions, times, extended file system attributes). If either the source
    file does not exist *or* the target file already exists, an exception will
    be raised.

    Such regular expression may be either a string *or* `Pattern` (i.e.,
    compiled regular expression object), while such substitution may be either a
    string *or* function. This function accepts the same optional keyword
    arguments as `re.sub()`.
    '''
    assert isinstance(filename_source, str),\
        '"{}" not a string.'.format(filename_source)
    assert isinstance(filename_target, str),\
        '"{}" not a string.'.format(filename_target)
    assert len(filename_source), 'Source filename empty.'
    assert len(filename_target), 'Target filename empty.'

    # Log such copy.
    loggers.log_info(
        'Munging file "%s" to "%s".', filename_source, filename_target)

    # Raise an exception unless the source file exists.
    die_unless_file(filename_source)

    # Raise an exception if the target file already exists.
    die_if_file(filename_target)

    # For atomicity, incrementally write to a temporary file rather than the
    # desired target file *BEFORE* moving the former to the latter. This
    # obscures such writes from other threads and/or processes, avoiding
    # potential race conditions elsewhere.
    with open_for_text_writing_temporary() as file_target_temp:
        with open_for_text_reading(filename_source) as file_source:
            for line in file_source:
                file_target_temp.write(re.sub(regex, substitution, line))

        # Copy all metadata (e.g., permissions) from such source to target file
        # *BEFORE* moving the latter, avoiding potential race conditions and
        # security vulnerabilities elsewhere.
        shutil.copystat(filename_source, file_target_temp.name)

        # Move such temporary file to such target file.
        shutil.move(file_target_temp.name, filename_target)

# --------------------( WASTELANDS                         )--------------------
    # return path.isfile(filename)
    # Versus path.isfile()
    # ----------
    # This function fundamentally differs from the stock `path.isfile()` function.
    # Whereas the latter returns True only for non-special files and hence False
    # for all non-directory special files (e.g., device nodes, symbolic links),
    # this function returns True for for *all* non-directory files regardless of
    # whether such files are special or not.
    #
    # Why? Because this function complies with POSIX semantics, whereas
    # `path.isfile()` does *not*. Under POSIX, it is largely irrelevant whether a non-directory
    # file is special or not; it simply matters whether such file is a directory
    # or not. For example, the external command `rm` removes only non-directory
    # files and the external command `rmdir` removes only empty directories.
    # '''
    # assert isinstance(filename, str), '"{}" not a string.'.format(filename)
    # return path.exists(filename) and not path.isdir(filename)

    # Such file will be copied in a manner preserving some but *not* all metadata,
    # in accordance with standard POSIX behaviour. Specifically, the permissions
    # but *not* owner, group, or times of such file

    # Such file will be copied in a manner maximally preserving metadata (e.g.,
    # owner, group, permissions, times, extended file system attributes
    # If such source file is a symbolic link, such link rather than the target
    # file of such link will be copied.
    # If such file does *not* exist, an exception is raised.
