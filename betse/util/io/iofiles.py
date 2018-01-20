#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level non-directory file content facilities.

See Also
----------
:mod:`betse.util.path.files`
    Low-level non-directory filename facilities.
'''

# ....................{ IMPORTS                            }....................
import re, shutil
from betse.util.io.log import logs
from betse.util.type.types import type_check, FileType, SequenceTypes
from io import BufferedIOBase, TextIOWrapper

# ....................{ GLOBALS                            }....................
READLINE_EOF = ''
'''
String returned by the :class:`file.readline` method of all :class:`file`-like
objects signifying end-of-file (EOF).

By design, this method is guaranteed to unambiguously return this string if and
only if no additional lines exist to read. For readability, testing against this
human-readable global is strongly preferable to testing against the equally
valid empty string.

See Also
----------
https://docs.python.org/3/tutorial/inputoutput.html#methods-of-file-objects
    Further details.
'''

# ....................{ GETTERS                            }....................
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

    with reading_chars(filename=filename, encoding=encoding) as text_file:
        return text_file.read()

# ....................{ GETTERS ~ mode                     }....................
@type_check
def get_mode_write_chars(is_overwritable: bool = False) -> str:
    '''
    Mode string suitable for opening a file handle for character-oriented
    writing via the ``mode`` parameter to the :func:`open` builtin.

    This low-level I/O function is principally intended to be called by
    higher-level I/O functions (e.g., :func:`writing_chars`).

    Parameters
    ----------
    is_overwritable : optional[bool]
        ``True`` if overwriting this file when this file already exists *or*
        ``False`` if raising an exception when this file already exists.
        Defaults to ``False`` for safety.

    Returns
    ----------
    str
        If the ``is_overwritable`` parameter is:
        * ``True``, this is ``"xt"``, raising exceptions when attempting to
          write files that already exist with this mode.
        * ``False``, this is ``"wt"``, silently overwriting such files.
    '''

    return 'wt' if is_overwritable else 'xt'


@type_check
def get_mode_write_bytes(is_overwritable: bool = False) -> str:
    '''
    Mode string suitable for opening a file handle for byte-oriented writing via
    the ``mode`` parameter to the :func:`open` builtin.

    This low-level I/O function is principally intended to be called by
    higher-level I/O functions (e.g., :func:`writing_bytes`).

    Parameters
    ----------
    is_overwritable : optional[bool]
        ``True`` if overwriting this file when this file already exists *or*
        ``False`` if raising an exception when this file already exists.
        Defaults to ``False`` for safety.

    Returns
    ----------
    str
        If the ``is_overwritable`` parameter is:
        * ``True``, this is ``"xb"``, raising exceptions when attempting to
          write files that already exist with this mode.
        * ``False``, this is ``"wb"``, silently overwriting such files.
    '''

    return 'wb' if is_overwritable else 'xb'

# ....................{ READERS ~ context                  }....................
@type_check
def reading_bytes(filename: str) -> BufferedIOBase:
    '''
    Open and return a filehandle suitable for reading the binary file with the
    passed filename, transparently decompressing this file if the filetype of
    this filename is that of a supported archive format.

    This function returns a :class:`file`-like object suitable for use wherever
    the :func:`open` builtin is callable (e.g., in ``with`` statements).

    Parameters
    ----------
    filename : str
        Relative or absolute path of the binary file to be read. If this
        filename is suffixed by a supported archive filetype (i.e., if the
        :func:`betse.util.path.archives.is_filetype` function returns ``True``
        for this filename), the returned filehandle automatically reads the
        decompressed rather than compressed byte contents of this file.

    Returns
    ----------
    BufferedIOBase
        :class:`file`-like object encapsulating this opened file.

    Raises
    ----------
    BetseFileException
        If this filename is *not* that of an existing file.
    '''

    # Avoid circular import dependencies.
    from betse.util.path import archives, files

    # Log this I/O operation.
    logs.log_debug('Reading bytes: %s', filename)

    # Raise an exception unless this file exists.
    files.die_unless_file(filename)

    # If this file is compressed, open and return a file handle reading
    # decompressed bytes from this file.
    if archives.is_filetype(filename):
        return archives.read_bytes(filename)
    # Else, this file is uncompressed. Open and return a typical file handle
    # reading bytes from this file.
    else:
        return open(filename, mode='rb')


@type_check
def reading_chars(filename: str, encoding: str = 'utf-8') -> TextIOWrapper:
    '''
    Open and return a filehandle suitable for reading the plaintext file with
    the passed filename encoded with the passed encoding.

    This function returns a :class:`file`-like object suitable for use wherever
    the :func:`open` builtin is callable (e.g., in ``with`` statements).

    Parameters
    ----------
    filename : str
        Relative or absolute path of the plaintext text to be read.
    encoding : optional[str]
        Name of the encoding to be used. Defaults to UTF-8.

    Returns
    ----------
    TextIOWrapper
        :class:`file`-like object encapsulating the opened file.
    '''

    # Avoid circular import dependencies.
    from betse.util.path import files

    # Log this I/O operation.
    logs.log_debug('Reading chars: %s', filename)

    # Raise an exception unless this file exists.
    files.die_unless_file(filename)

    # Open this file.
    return open(filename, mode='rt', encoding=encoding)

# ....................{ WRITERS                            }....................
@type_check
def write_str_to_filename(text: str, output_filename: str, **kwargs) -> None:
    '''
    Overwrite the plaintext output file with the passed filename with the
    passed string.

    Parameters
    ----------
    text : FileType
        Contents to be written to this plaintext file.
    output_filename : str
        Absolute or relative path of the plaintext file to be overwritten.

    All remaining keyword arguments are passed as is to the
    :func:`writing_chars` context manager.
    '''

    with writing_chars(filename=output_filename, **kwargs) as output_file:
        output_file.write(text)


@type_check
def write_file_to_filename(
    input_file: FileType, output_filename: str, **kwargs) -> None:
    '''
    Overwrite the output file with the passed filename with the contents of the
    passed open readable file-like object.

    To copy the entirety of this input file, the seek position of this file-like
    object is reset to the first byte of this file. Callers preferring that this
    position be preserved should do so manually.

    Parameters
    ----------
    input_file : FileType
        Open readable file-like object to copy.
    output_filename : str
        Absolute or relative path of the file to copy into.

    All remaining keyword arguments are passed as is to the
    :func:`writing_chars` context manager.
    '''

    #FIXME: Actually, wouldn't it be preferable to conditionally:
    #
    #* Call the writing_bytes() function if this "input_file" is currently open
    #  in a byte-orintied manner.
    #* Call the writing_chars() function in all other cases as a fallback.
    #
    #How, exactly, would one perform the former check?

    # For portability, open this output file for character-oriented writing.
    # Opening this file for byte-oriented writing would technically be more
    # efficient but fail for all character-oriented file-like objects (e.g.,
    # instances of the standard "io.StringIO" class).
    with writing_chars(filename=output_filename, **kwargs) as output_file:
        write_file_to_file(input_file=input_file, output_file=output_file)


@type_check
def write_file_to_file(input_file: FileType, output_file: FileType) -> None:
    '''
    Overwrite the passed open writable file-like object with the contents of the
    passed open readable file-like object.

    To copy the entirety of this input file, the seek positions of both this
    objects are reset to the first bytes of their corresponding files. Callers
    requiring that these positions be preserved should do so manually.

    Parameters
    ----------
    input_file : FileType
        Open readable file-like object to copy.
    output_file : FileType
        Open writable file-like object to copy into.
    '''

    # Reset the seek positions of both objects.
    input_file.seek(0)
    output_file.seek(0)

    # Safely copy in a manner avoiding space (i.e., memory) exhaustion.
    # print('types: {}, {}'.format(type(input_file), type(output_file)))
    shutil.copyfileobj(input_file, output_file)

# ....................{ WRITERS ~ context                  }....................
@type_check
def writing_bytes(
    filename: str, is_overwritable: bool = False) -> BufferedIOBase:
    '''
    Open and return a filehandle suitable for writing the binary file with the
    passed filename, transparently compressing this file if the filetype of
    this filename is that of a supported archive format.

    This function returns a file-like object suitable for use wherever the
    :func:`open` builtin is callable (e.g., in ``with`` statements).

    Parameters
    ----------
    filename : str
        Relative or absolute path of the binary file to be written. If this
        filename is suffixed by a supported archive filetype (i.e., if the
        :func:`betse.util.path.archives.is_filetype` function returns ``True``
        for this filename), the returned filehandle automatically writes the
        compressed rather than uncompressed byte contents of this file.
    is_overwritable : optional[bool]
        ``True`` if overwriting this file when this file already exists *or*
        ``False`` if raising an exception when this file already exists.
        Defaults to ``False`` for safety.

    Returns
    ----------
    BufferedIOBase
        File-like object encapsulating this opened file.
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


@type_check
def writing_chars(
    filename: str,
    is_overwritable: bool = False,
    encoding: str = 'utf-8',
) -> TextIOWrapper:
    '''
    Open and return a filehandle suitable for writing the plaintext file with
    the passed filename encoded with the passed encoding.

    This function returns a file-like object suitable for use wherever the
    :func:`open` builtin is callable (e.g., in ``with`` statements).

    Parameters
    ----------
    filename : str
        Relative or absolute path of the plaintext file to be written.
    is_overwritable : optional[bool]
        ``True`` if overwriting this file when this file already exists *or*
        ``False`` if raising an exception when this file already exists.
        Defaults to ``False`` for safety.
    encoding : optional[str]
        Name of the encoding to be used. Defaults to UTF-8.

    Returns
    ----------
    TextIOWrapper
        File-like object encapsulating this opened file.
    '''

    # Avoid circular import dependencies.
    from betse.util.path import dirs, paths

    # Log this I/O operation.
    logs.log_debug('Writing chars: %s', filename)

    # If this file is *NOT* overwritable, raise an exception if this path
    # already exists.
    if not is_overwritable:
        paths.die_if_path(filename)

    # Create the parent directory of this file if needed.
    dirs.make_parent_unless_dir(filename)

    # Mode with which to open this file for character-oriented writing.
    mode = get_mode_write_chars(is_overwritable)

    # Open and return a file handle writing uncompressed bytes to this file.
    return open(filename, mode=mode, encoding=encoding)

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

    This function implements the equivalent of the ``sed`` line processor in a
    pure-Python manner requiring no external commands or additional
    dependencies. This is a good thing.

    This source file will *not* be changed. This target file will be written in
    an atomic manner maximally preserving source metadata (e.g., owner, group,
    permissions, times, extended file system attributes). If either the source
    file does not exist *or* the target file already exists, an exception is
    raised.

    These regular expressions may be either strings *or* instances of the
    :class:`Pattern` class (i.e., compiled regular expression object), while
    these substitutions may be either strings *or* functions. See the
    "Arguments" section below for how this function expects these objects to be
    passed.

    This function accepts the same optional keyword arguments as the standard
    :func:`re.sub` function.

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
    from betse.util.path import files, temps

    # Log this substitution.
    if filename_source == filename_target:
        logs.log_debug(
            'Munging file "%s" in-place.', filename_source)
    else:
        logs.log_debug(
            'Munging file "%s" to "%s".', filename_source, filename_target)

    # Raise an exception unless the source file exists.
    files.die_unless_file(filename_source)

    # Raise an exception if the target file already exists.
    files.die_if_file(filename_target)

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
        with reading_chars(filename_source) as file_source:
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
