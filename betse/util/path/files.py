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
from betse.util.type import types
from os import path
import os, re, shutil, tempfile

# ....................{ EXCEPTIONS ~ unless                }....................
def die_unless_file(pathname: str) -> None:
    '''
    Raise an exception unless the passed non-directory file exists _after_
    following symbolic links.
    '''
    if not is_file(pathname):
        raise BetseExceptionFile(
            'File "{}" not found or unreadable.'.format(pathname))

# ....................{ EXCEPTIONS ~ if                    }....................
def die_if_file(pathname: str) -> None:
    '''
    Raise an exception if the passed non-directory file exists _after_ following
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

    **Why?** Because this function complies with POSIX semantics, whereas
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
    assert types.is_str_nonempty(pathname),\
        types.assert_is_nonstr_nonempty(pathname, 'pathname')
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
    assert types.is_str_nonempty(filename_source),\
        types.assert_is_nonstr_nonempty(filename_source, 'source filename')
    assert types.is_str_nonempty(filename_target),\
        types.assert_is_nonstr_nonempty(filename_target, 'target filename')

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
    assert types.is_str_nonempty(filename),\
        types.assert_is_nonstr_nonempty(filename, 'filename')

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
    assert types.is_str_nonempty(filename),\
        types.assert_is_nonstr_nonempty(filename, 'filename')

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
    assert types.is_str_nonempty(filename),\
        types.assert_is_nonstr_nonempty(filename, 'filename')

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

def substitute_substrings(
    filename_source: str,
    filename_target: str,
    substitutions,
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
    substitutions : sequence_nonstring
        Non-string sequence (e.g., list, tuple) of non-string sequences of
        length 2 (i.e., pairs), whose first element is a regular expression and
        whose second element is the substitution to be performed for all
        substrings in the source file matching that regular expression.
    '''
    assert types.is_str_nonempty(filename_source),\
        types.assert_is_nonstr_nonempty(filename_source, 'source filename')
    assert types.is_str_nonempty(filename_target),\
        types.assert_is_nonstr_nonempty(filename_target, 'target filename')
    assert types.is_sequence_nonstr_nonempty(substitutions),\
        types.assert_is_not_sequence_nonstr_nonempty(
            substitutions, 'regular expression substitution pairs')

    # Log such substitution.
    if filename_source == filename_target:
        loggers.log_info(
            'Munging file "%s" in-place.', filename_source)
    else:
        loggers.log_info(
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

    #FIXME: Such functionality is probably quite useful, where at least one
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

# --------------------( WASTELANDS                         )--------------------
#FUXME: Correct documentation here and below to reflect the passing of multiple
#regexes and substitutions.

                    # print('Yay!')
                    # print('Line: {}'.format(line))
                    # if re.search(regex, line, **kwargs) is not None:
                    #     is_line_matches = True
                    #     loggers.log_info('Line "%s" matches!', line)
                    # if re.search(
                    #     # r'^(\s*turn all plots off)',
                    #     r'turn all plots off', line) is not None:
                    #     is_line_matches = True
                    #     loggers.log_info('Line "%s" really matches!', line)
                    # line = re.sub(regex, substitution, line, **kwargs)

        # (regex, substitution)
        # for substitution_pair in substitutions
        # for regex, substitution in substitution_pair

#     assert types.is_sequence_nonstr_nonempty(regexes),\
#         types.assert_is_not_sequence_nonstr_nonempty(
#             regexes, 'regular expressions')
#     assert types.is_sequence_nonstr_nonempty(substitutions),\
#         types.assert_is_not_sequence_nonstr_nonempty(
#             substitutions, 'substitutions')
#     assert len(regexes) == len(substitutions),\
#         '{} regular expressions unequal to {} substitutions.'.format(
#             len(regexes), len(substitutions))
# #The number of passed regular expressions _must_ equal the number of passed substitutions.
#     # return path.isfile(filename)
#     # Versus path.isfile()
#     # ----------
#     # This function fundamentally differs from the stock `path.isfile()` function.
#     # Whereas the latter returns True only for non-special files and hence False
#     # for all non-directory special files (e.g., device nodes, symbolic links),
#     # this function returns True for for *all* non-directory files regardless of
#     # whether such files are special or not.
#     #
#     # Why? Because this function complies with POSIX semantics, whereas
#     # `path.isfile()` does *not*. Under POSIX, it is largely irrelevant whether a non-directory
#     # file is special or not; it simply matters whether such file is a directory
#     # or not. For example, the external command `rm` removes only non-directory
#     # files and the external command `rmdir` removes only empty directories.
#     # '''
#     # assert isinstance(filename, str), '"{}" not a string.'.format(filename)
#     # return path.exists(filename) and not path.isdir(filename)

    # Such file will be copied in a manner preserving some but *not* all metadata,
    # in accordance with standard POSIX behaviour. Specifically, the permissions
    # but *not* owner, group, or times of such file

    # Such file will be copied in a manner maximally preserving metadata (e.g.,
    # owner, group, permissions, times, extended file system attributes
    # If such source file is a symbolic link, such link rather than the target
    # file of such link will be copied.
    # If such file does *not* exist, an exception is raised.
