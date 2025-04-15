#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2025 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level directory facilities.
'''

#FIXME: Splint into multiple submodules for maintainability by defining a new
#"betse.util.path.directory" subpackage to store these submodules. Note that
#this subpackage is intentionally *NOT* named "dir", as doing so would conflict
#with the standard dir() builtin.

# ....................{ IMPORTS                            }....................
import os, shutil, time
from betse.exceptions import BetseDirException, BetsePathException
from betse.util.io.log import logs
from betse.util.type.enums import make_enum
from betse.util.type.types import (
    type_check,
    CallableTypes,
    GeneratorType,
    IterableOrNoneTypes,
    NumericSimpleTypes,
    SequenceTypes,
)
from os import path as os_path

# ....................{ ENUMERATIONS                       }....................
DirOverwritePolicy = make_enum(
    class_name='DirOverwritePolicy',
    member_names=(
        'HALT_WITH_EXCEPTION',
        'SKIP_WITH_WARNING',
        'OVERWRITE',
    ),
)
'''
Enumeration of all supported types of **directory overwrite policies** (i.e.,
strategies for handling existing paths visited by recursive directory
operations).

Attributes
----------
HALT_WITH_EXCEPTION : enum
    Policy raising a fatal exception if any target path already exists. This
    constitutes the strictest and thus safest such policy.
SKIP_WITH_WARNING : enum
    Policy ignoring (i.e., skipping) each existing target path with a non-fatal
    warning. This policy strikes a comfortable balance between strictness and
    laxness and is thus the recommended default.
OVERWRITE : enum
    Policy silently overwriting each existing target path. This constitutes the
    laxest and thus riskiest such policy.
'''

# ....................{ GLOBALS                            }....................
# There exist only two possible directory separators for all modern platforms.
# Hence, this reliably suffices with no error handling required.
SEPARATOR_REGEX = r'/' if os_path.sep == '/' else r'\\'
'''
Regular expression matching the directory separator specific to the current
platform.

Specifically, under:

* Microsoft Windows, this is `\\\\`.
* All other platforms, this is `/`.
'''

# ....................{ EXCEPTIONS                         }....................
@type_check
def die_if_dir(*dirnames: str) -> None:
    '''
    Raise an exception if any of the directories with the passed dirnames
    exist.

    Parameters
    ----------
    dirnames: Tuple[str]
        Tuple of the absolute or relative dirnames of all directories to be
        validated.

    Raises
    ----------
    BetseDirException
        If any passed directory exists.
    '''

    for dirname in dirnames:
        if is_dir(dirname):
            raise BetseDirException(
                'Directory "{}" already exists.'.format(dirname))


@type_check
def die_unless_dir(*dirnames: str) -> None:
    '''
    Raise an exception unless all of the directories with the passed dirnames
    exist.

    Equivalently, this function raises an exception if any of these directories
    do *not* exist.

    Parameters
    ----------
    dirnames: Tuple[str]
        Tuple of the absolute or relative dirnames of all directories to be
        validated.

    Raises
    ----------
    BetseDirException
        If any passed directory does *not* exist.
    '''

    for dirname in dirnames:
        if not is_dir(dirname):
            raise BetseDirException(
                'Directory "{}" not found or unreadable.'.format(dirname))


@type_check
def dir_or_die(dirname: str) -> str:
    '''
    Passed dirname if a directory with this dirname exists *or* raise an
    exception otherwise (i.e., if this directory does *not* exist).

    Parameters
    ----------
    dirname : str
        Absolute or relative dirname of the directory to be validated.

    Returns
    ----------
    str
        Passed dirname if this directory exists.

    Raises
    ----------
    BetseDirException
        If this directory does *not* exist.
    '''

    # If this directory is not found, fail.
    die_unless_dir(dirname)

    # Else, this directory exists. Return this directory's pathname.
    return dirname


@type_check
def join_or_die(*pathnames: str) -> str:
    '''
    Pathname of the directory produced by joining (i.e., concatenating) the
    passed pathnames with the directory separator specific to the current
    platform if this directory exists *or* raise an exception otherwise (i.e.,
    if this directory does *not* exist).

    This higher-level function is a convenience wrapper encapsulating both the
    lower-level :func:`betse.util.path.pathnames.join` and
    :func:`die_unless_dir` functions.
    '''

    # Avoid circular import dependencies.
    from betse.util.path.pathnames import join

    # Dirname producing by joining these pathnames.
    dirname = join(*pathnames)

    # If this directory does *NOT* exist, raise an exception.
    die_unless_dir(dirname)

    # Return this dirname.
    return dirname

# ....................{ EXCEPTIONS ~ subdir                }....................
@type_check
def die_if_subdir(parent_dirname: str, child_dirname: str) -> None:
    '''
    Raise an exception if the child directory with the passed dirname actually
    is a child (i.e., subdirectory) of the parent directory with the passed
    dirname.

    Parameters
    ----------
    parent_dirname: str
        Absolute or relative dirname of the parent directory to be validated.
    child_dirname: str
        Absolute or relative dirname of the child directory to be validated.

    Raises
    ----------
    BetseDirException
        If this child directory is *not* actually a child (i.e., subdirectory)
        of this parent directory.
    '''

    if is_subdir(parent_dirname, child_dirname):
        raise BetseDirException(
            '"{}" is a subdirectory of "{}".'.format(
                child_dirname, parent_dirname))

# ....................{ EXCEPTIONS ~ parent                }....................
@type_check
def die_unless_parent_dir(pathname: str) -> None:
    '''
    Raise an exception unless the parent directory of the passed path exists.
    '''

    # Avoid circular import dependencies.
    from betse.util.path import pathnames

    # If the parent directory of this path does *NOT* exist, raise an
    # exception.
    die_unless_dir(pathnames.get_dirname(pathname))

# ....................{ TESTERS                            }....................
@type_check
def is_dir(dirname: str) -> bool:
    '''
    ``True`` only if the directory with the passed dirname exists.

    Parameters
    -----------
    dirname : str
        Absolute or relative dirname of the directory to be tested.

    Returns
    -----------
    bool
        ``True`` only if this directory exists.
    '''

    # One-liners for terrible vengeance.
    return os_path.isdir(dirname)


@type_check
def is_empty(dirname: str) -> bool:
    '''
    ``True`` only if the directory with the passed dirname is **empty** (i.e.,
    contains no files or subdirectories).

    Parameters
    -----------
    dirname : str
        Absolute or relative dirname of the directory to be tested.

    Returns
    -----------
    bool
        ``True`` only if this directory is empty.

    Raises
    -----------
    BetseDirException
        If the passed dirname either is not an existing path *or* is but is not
        a directory.
    '''

    # If this path is *NOT* a directory, raise an exception.
    die_unless_dir(dirname)

    # Return true only if this directory contains no child paths.
    return not os.listdir(dirname)


@type_check
def is_subdir(parent_dirname: str, child_dirname: str) -> bool:
    '''
    ``True`` only if the child directory with the passed dirname actually is a
    child (i.e., subdirectory) of the parent directory with the passed dirname.

    Parameters
    -----------
    parent_dirname: str
        Absolute or relative dirname of the parent directory to be tested.
    child_dirname: str
        Absolute or relative dirname of the child directory to be tested.

    Returns
    -----------
    bool
        ``True`` only if this child directory actually is a child (i.e.,
        subdirectory) of this parent directory.

    See Also
    -----------
    https://stackoverflow.com/a/37095733/2809027
        StackOverflow answer strongly inspiring this implementation.
    '''

    # Avoid circular import dependencies.
    from betse.util.os.brand.windows import is_windows
    from betse.util.path import pathnames

    # Canonicalized child and parent dirnames, conditionally resolving both
    # relative dirnames and symbolic links as needed.
    #
    # Note that the subsequently called os_path.commonpath() function
    # explicitly requires these dirnames to be canonicalized, raising a
    # "ValueError" exception if this is *NOT* the case.
    parent_dirname = pathnames.canonicalize(parent_dirname)
    child_dirname  = pathnames.canonicalize(child_dirname)

    # Longest common dirname shared between these dirnames if any *OR* the root
    # directory otherwise (e.g., "/" under Linux), defaulting to the empty
    # string to ensure the comparison below fails if this is *NOT* set below.
    common_dirname = ''

    # Attempt to decide this longest common dirname.
    try:
        common_dirname = os_path.commonpath((parent_dirname, child_dirname))
    # If commonpath() raised a "ValueError" exception *AND*...
    except ValueError as exception:
        #FIXME: Consider raising an upstream issue against the CPython issue
        #tracker. This behaviour is unexpected and hostile to platform-portable
        #resilience throughout the Python ecosystem.

        # The current platform is Microsoft Windows *AND* this exception message
        # is a well-known message raised by the Windows-specific standard
        # "ntpath" module indicating the two passed pathnames to reside on
        # different drives (e.g., "C:/" and "D:/"), silently ignore this
        # exception. This Windows-specific behaviour is fundamentally broken.
        # Since two pathnames residing on different drives have no common path,
        # the Windows-specific commonpath() implementation should return false
        # rather than raise an exception. Force this to be the case here!
        if is_windows() and str(exception) == "Paths don't have the same drive":
            pass
        # Else, re-raise this unexpected exception.
        else:
            raise
    # Else, the commonpath() builtin raises a "ValueError" exception *AND*...

    # Return true only if this dirname is this parent's canonicalized dirname.
    return parent_dirname == common_dirname

# ....................{ GETTERS                            }....................
@type_check
def get_parent_dir_last(pathname: str) -> str:
    '''
    Absolute pathname of the most deeply nested existing parent directory of
    the path with the passed pathname.

    Since the passed pathname is required to be absolute *and* since the root
    directory (e.g., ``/`` on POSIX-compatible platforms) always exists, this
    function is guaranteed to return the pathname of an existing directory.

    Parameters
    -----------
    pathname : str
        Absolute pathname of the path to inspect.

    Returns
    -----------
    str
        Absolute pathname of the most deeply nested existing parent directory
        of this path.

    Examples
    -----------
    Given an existing directory ``/the/garden/of/earthly/`` and non-existing
    subdirectory ``/the/garden/of/earthly/delights/`` of that directory:

        >>> from betse.util.path import dirs
        >>> dirs.get_parent_dir_last('/the/garden/of/earthly/delights/')
        '/the/garden/of/earthly'
    '''

    # Avoid circular import dependencies.
    from betse.util.path import pathnames

    # If this pathname is relative rather than absolute, raise an exception.
    pathnames.die_if_relative(pathname)

    # Current directory component of this path to inspect the existence of.
    dirname_component = pathname

    # While this directory component is *NOT* an existing directory...
    while not is_dir(dirname_component):
        # Reduce this directory component to its parent directory.
        dirname_component = pathnames.get_dirname(dirname_component)

    # Return this existing directory component.
    return dirname_component

# ....................{ GETTERS ~ mtime : recursive        }....................
#FIXME: Consider contributing the get_mtime_recursive_newest() function as an
#answer to the following StackOverflow question:
#    https://stackoverflow.com/questions/26498285/find-last-update-time-of-a-directory-includes-all-levels-of-sub-folders

# If the current platform provides the os.fwalk() function *AND* the os.stat()
# function implemented by this platform accepts directory handles (which is
# currently equivalent to testing if this platform is POSIX-compatible),
# optimize this recursion by leveraging these handles.
if hasattr(os, 'fwalk') and os.stat in os.supports_dir_fd:
    @type_check
    def get_mtime_recursive_newest(dirname: str) -> NumericSimpleTypes:

        # Log this recursion.
        logs.log_debug(
            'Recursively computing newest fwalk()-based mtime of: %s', dirname)

        # Return the maximum of all mtimes in the...
        return max(
            # Generator expression recursively yielding the maximum mtime of
            # all files and subdirectories of each subdirectory of this
            # directory including this directory.
            (
                # Maximum of this subdirectory's mtime and the mtimes of all
                # files in this subdirectory, whose filenames are produced by a
                # tuple comprehension joining this subdirectory's dirname to
                # each file's basename. By recursion, the mtimes of all
                # subdirectories of this subdirectory are implicitly computed
                # and hence need *NOT* be explicitly included here.
                #
                # For parity with the unoptimized get_mtime_recursive_newest()
                # implementation defined below as well to avoid unwanted
                # complications, symbolic links are *NOT* followed. Hence,
                # os.lstat() rather than os.stat() is intentionally called.
                max(
                    (os.fstat(parent_dir_fd).st_mtime,) + tuple(
                        os.lstat(
                            child_file_basename, dir_fd=parent_dir_fd).st_mtime
                        for child_file_basename in child_file_basenames
                    ),
                )
                # For the currently visited directory's dirname, a sequence of
                # the dirnames of all subdirectories of this directory, a
                # sequence of the filenames of all files in this directory, and
                # a directory handle to this directory such that errors emitted
                # by low-level functions called by os.walk() (e.g.,
                # os.listdir()) are *NOT* silently ignored...
                for _, _, child_file_basenames, parent_dir_fd in _fwalk(dirname)
            ),
        )

# Else, fallback to unoptimized recursion leveraging dirnames.
else:
    @type_check
    def get_mtime_recursive_newest(dirname: str) -> NumericSimpleTypes:

        # Log this recursion.
        logs.log_debug(
            'Recursively computing newest walk()-based mtime of: %s', dirname)

        # Return the maximum of all mtimes in the...
        return max(
            # Generator expression recursively yielding the maximum mtime of
            # all files and subdirectories of each subdirectory of this
            # directory including this directory.
            (
                # Maximum of this subdirectory's mtime and the mtimes of all
                # files in this subdirectory, whose filenames are produced by a
                # tuple comprehension joining this subdirectory's dirname to
                # each file's basename. By recursion, the mtimes of all
                # subdirectories of this subdirectory are implicitly computed
                # and hence need *NOT* be explicitly included here.
                max(
                    (os_path.getmtime(parent_dirname),) + tuple(
                        os_path.getmtime(os_path.join(
                            parent_dirname, child_file_basename))
                        for child_file_basename in child_file_basenames
                    ),
                )
                # For the currently visited directory's dirname, a sequence of
                # the dirnames of all subdirectories of this directory, and a
                # sequence of the filenames of all files in this directory such
                # that errors emitted by low-level functions called by
                # os.walk() (e.g., os.listdir()) are *NOT* silently ignored...
                for parent_dirname, _, child_file_basenames in _walk(dirname)
            ),
        )

# Docstring dynamically set for the getter defined above.
get_mtime_recursive_newest.__doc__ = '''
    Most recent **recursive mtime** (i.e., recursively calculated modification
    time) in seconds of all paths in the set of this directory and all files and
    subdirectories transitively reachable from this directory *witthout*
    following symbolic links to directories.

    Caveats
    -----------
    Note that symbolic links whose transitive targets are:

    * Directories are *not* followed, avoiding infinite recursion in edge cases
      (e.g., directories containing symbolic links to themselves). These links
      are effectively non-directory files. Since the mtime of a symbolic link
      is typically its ctime (i.e., creation time), such links thus contribute
      their own ctime to this calculation.
    * Files are followed, improving the accuracy of this calculation.

    Parameters
    -----------
    dirname : str
        Absolute or relative path of the directory to inspect.

    Returns
    -------
    NumericSimpleTypes
        Most recent mtime in seconds of this directory calculated recursively
        as either an integer or float, depending on the boolean returned by the
        platform-specific :func:`os.stat_float_times` function.

    Raises
    ------
    OSError
        If either this directory *or* any file or subdirectory reachable from
        this directory is unreadable by the current user (e.g., due to
        insufficient permissions).
    '''

# ....................{ COPIERS                            }....................
@type_check
def copy_dir_into_dir(
    src_dirname: str, trg_dirname: str, *args, **kwargs) -> None:
    '''
    Recursively copy the source directory with the passed dirname to a
    subdirectory of the target directory with the passed dirname whose basename
    is that of this source directory.

    Parameters
    -----------
    src_dirname : str
        Absolute or relative dirname of the source directory to be copied from.
    trg_dirname : str
        Absolute or relative dirname of the target directory to be copied into.

    All remaining parameters are passed as is to the :func:`copy` function.

    See Also
    ----------
    :func:`copy_dir`
        Further details.

    Examples
    ----------
        >>> from betse.util.path import dirs
        >>> dirs.copy_into_dir('/usr/src/linux/', '/tmp/')
        >>> dirs.is_dir('/tmp/linux/')
        True
    '''

    # Avoid circular import dependencies.
    from betse.util.path import pathnames

    # Basename of this source directory.
    src_basename = pathnames.get_basename(src_dirname)

    # Absolute or relative dirname of the target directory to copy to.
    trg_subdirname = pathnames.join(trg_dirname, src_basename)

    # Recursively copy this source to target directory.
    copy_dir(
        *args, src_dirname=src_dirname, trg_dirname=trg_subdirname, **kwargs)


@type_check
def copy_dir(
    # Mandatory parameters.
    src_dirname: str,
    trg_dirname: str,

    # Optional parameters.
    overwrite_policy: DirOverwritePolicy = (
        DirOverwritePolicy.HALT_WITH_EXCEPTION),
    ignore_basename_globs: IterableOrNoneTypes = None,
) -> None:
    '''
    Recursively copy the source directory with the passed dirname into the
    target directory with the passed dirname.

    For generality:

    * All nonexistent parents of the target directory will be recursively
      created, mimicking the action of the ``mkdir -p`` shell command on
      POSIX-compatible platforms in a platform-agnostic manner.
    * All symbolic links in the source directory will be preserved (i.e.,
      copied as is rather than their transitive targets copied instead).

    Parameters
    -----------
    src_dirname : str
        Absolute or relative dirname of the source directory to be copied from.
    trg_dirname : str
        Absolute or relative dirname of the target directory to be copied to.
    overwrite_policy : DirOverwritePolicy
        **Directory overwrite policy** (i.e., strategy for handling existing
        paths to be overwritten by this copy). Defaults to
        :attr:`DirOverwritePolicy.HALT_WITH_EXCEPTION`, raising an exception if
        any target path already exists.
    ignore_basename_globs : optional[IterableTypes]
        Iterable of shell-style globs (e.g., ``('*.tmp', '.keep')``) matching
        the basenames of all paths transitively owned by this source directory
        to be ignored during recursion and hence neither copied nor visited.
        If non-``None`` and the ``overwrite_policy`` parameter is
        :attr:`DirOverwritePolicy.OVERWRITE`, this iterable is ignored and a
        non-fatal warning is logged. Defaults to ``None``, in which case *all*
        paths transitively owned by this source directory are unconditionally
        copied and visited.

    Raises
    -----------
    BetseDirException
        If either:

        * The source directory does *not* exist.
        * The target directory is a subdirectory of the source directory.
          Permitting this edge case induces non-trivial issues, including
          infinite recursion from within the musty entrails of the
          :mod:`distutils` package (e.g., due to relative symbolic links).
        * The passed ``overwrite_policy`` parameter is
          :attr:`DirOverwritePolicy.HALT_WITH_EXCEPTION` *and* one or more
          subdirectories of the target directory already exist that are also
          subdirectories of the source directory. For safety, this function
          always preserves rather than overwrites existing target
          subdirectories.

    See Also
    -----------
    https://stackoverflow.com/a/22588775/2809027
        StackOverflow answer strongly inspiring this function's
        :attr:`DirOverwritePolicy.SKIP_WITH_WARNING` implementation.
    '''

    # Log this copy.
    logs.log_debug('Copying directory: %s -> %s', src_dirname, trg_dirname)

    # If the source directory does *NOT* exist, raise an exception.
    die_unless_dir(src_dirname)

    # If the target directory is a subdirectory of the source directory, raise
    # an exception. Permitting this edge case provokes issues, including
    # infinite recursion from within the musty entrails of the "distutils"
    # codebase (possibly due to relative symbolic links).
    die_if_subdir(parent_dirname=src_dirname, child_dirname=trg_dirname)

    # If passed an iterable of shell-style globs matching ignorable basenames,
    # convert this iterable into a predicate function of the form required by
    # the shutil.copytree() function. Specifically, this function accepts the
    # absolute or relative pathname of an arbitrary directory and an iterable
    # of the basenames of all subdirectories and files directly in this
    # directory; this function returns an iterable of the basenames of all
    # subdirectories and files in this directory to be ignored. This signature
    # resembles:
    #
    #     def ignore_basename_func(
    #         parent_dirname: str,
    #         child_basenames: IterableTypes) -> IterableTypes
    ignore_basename_func = None
    if ignore_basename_globs is not None:
        ignore_basename_func = shutil.ignore_patterns(*ignore_basename_globs)

    # If either:
    # * Raising a fatal exception if any target path already exists *OR*...
    # * Overwriting this target directory with this source directory...
    #
    # Then the standard shutil.copytree() function applies to this use case.
    if overwrite_policy in _DIR_OVERWRITE_POLICIES_COPYTREE:
        # Dictionary of all keyword arguments to pass to shutil.copytree(),
        # preserving symbolic links as is.
        copytree_kwargs = {
            'symlinks': True,
        }

        # If raising a fatal exception if any target path already exists, do so.
        # While we could defer to the exception raised by the shutil.copytree()
        # function for this case, this exception's message erroneously refers to
        # this directory as a file and is hence best avoided as unreadable:
        #     [Errno 17] File exists: 'sample_sim'  # <-- lolbro! useless.
        if overwrite_policy is DirOverwritePolicy.HALT_WITH_EXCEPTION:
            die_if_dir(trg_dirname)
        # Else, we are overwriting this target directory with this source
        # directory. In this case, silently accept this target directory if this
        # directory already exists.
        else:
            copytree_kwargs['dirs_exist_ok'] = True

        # If ignoring basenames, inform shutil.copytree() of these basenames.
        if ignore_basename_func is not None:
            copytree_kwargs['ignore'] = ignore_basename_func

        # Recursively copy this source to target directory. To avoid silently
        # overwriting all conflicting target paths, the shutil.copytree()
        # rather than dir_util.copy_tree() function is called.
        shutil.copytree(src=src_dirname, dst=trg_dirname, **copytree_kwargs)

    #FIXME: Given how awesomely flexible the manual approach implemented below
    #is, we should probably consider simply rewriting the above two approaches
    #to reuse the exact same logic. It works. It's preferable. Let's reuse it.
    #FIXME: Actually, this is increasingly critical. Third-party functions
    #called above -- notably, the dir_util.copy_tree() function -- appear to
    #suffer critical edge cases. This can be demonstrated via the BETSEE GUI by
    #attempting to save an opened simulation configuration to a subdirectory of
    #itself, which appears to provoke infinite recursion from within the musty
    #depths of the "distutils" codebase. Of course, the implementation below
    #could conceivably suffer similar issues. If this is the case, this
    #function should explicitly detect attempts to recursively copy a source
    #directory into a subdirectory of itself and raise an exception.
    #FIXME: See the above FIXME comment addressing the infinite recursion issue
    #discussed here.

    # Else if logging a warning for each target path that already exists, do so
    # by manually implementing recursive directory copying. Sadly, Python
    # provides no means of doing so "out of the box."
    elif overwrite_policy is DirOverwritePolicy.SKIP_WITH_WARNING:
        # Avoid circular import dependencies.
        from betse.util.path import files, paths, pathnames
        from betse.util.type.iterable import sequences

        # Passed parameters renamed for disambiguity.
        src_root_dirname = src_dirname
        trg_root_dirname = trg_dirname

        # Basename of the top-level target directory to be copied to.
        trg_root_basename = pathnames.get_basename(src_root_dirname)

        # For the absolute pathname of each recursively visited source
        # directory, an iterable of the basenames of all subdirectories of this
        # directory, and an iterable of the basenames of all files of this
        # directory...
        for src_parent_dirname, subdir_basenames, file_basenames in _walk(
            src_root_dirname):
            # Relative pathname of the currently visited source directory
            # relative to the absolute pathname of this directory.
            parent_dirname_relative = pathnames.relativize(
                src_dirname=src_root_dirname, trg_pathname=src_parent_dirname)

            # If ignoring basenames...
            if ignore_basename_func is not None:
                # Sets of the basenames of all ignorable subdirectories and
                # files of this source directory.
                subdir_basenames_ignored = ignore_basename_func(
                    src_parent_dirname, subdir_basenames)
                file_basenames_ignored = ignore_basename_func(
                    src_parent_dirname, file_basenames)

                # If ignoring one or more subdirectories...
                if subdir_basenames_ignored:
                    # Log the basenames of these subdirectories.
                    logs.log_debug(
                        'Ignoring source "%s/%s" subdirectories: %r',
                        trg_root_basename,
                        parent_dirname_relative,
                        subdir_basenames_ignored)

                    # Remove these subdirectories from the original iterable.
                    # Since the os.walk() function supports in-place changes to
                    # this iterable, this iterable is modified via this less
                    # efficient function rather than efficient alternatives
                    # (e.g., set subtraction).
                    sequences.remove_items(
                        sequence=subdir_basenames,
                        items=subdir_basenames_ignored)

                # If ignoring one or more files...
                if file_basenames_ignored:
                    # Log the basenames of these files.
                    logs.log_debug(
                        'Ignoring source "%s/%s" files: %r',
                        trg_root_basename,
                        parent_dirname_relative,
                        file_basenames_ignored)

                    # Remove these files from the original iterable. Unlike
                    # above, we could technically modify this iterable via
                    # set subtraction: e.g.,
                    #
                    #     subdir_basenames -= subdir_basenames_ignored
                    #
                    # For orthogonality, preserve the above approach instead.
                    sequences.remove_items(
                        sequence=file_basenames,
                        items=file_basenames_ignored)

            # Absolute pathname of the corresponding target directory.
            trg_parent_dirname = pathnames.join(
                trg_root_dirname, parent_dirname_relative)

            # Create this target directory if needed.
            make_unless_dir(trg_parent_dirname)

            # For the basename of each non-ignorable file of this source
            # directory...
            for file_basename in file_basenames:
                # Absolute filenames of this source and target file.
                src_filename = pathnames.join(
                    src_parent_dirname, file_basename)
                trg_filename = pathnames.join(
                    trg_parent_dirname, file_basename)

                # If this target file already exists...
                if paths.is_path(trg_filename):
                    # Relative filename of this file. The absolute filename of
                    # this source or target file could be logged instead, but
                    # this relative filename is significantly more terse.
                    filename_relative = pathnames.join(
                        trg_root_basename,
                        parent_dirname_relative,
                        file_basename)

                    # Warn of this file being ignored.
                    logs.log_warning(
                        'Ignoring existing target file: %s', filename_relative)

                    # Ignore this file by continuing to the next.
                    continue

                # Copy this source to target file.
                files.copy(
                    src_filename=src_filename, trg_filename=trg_filename)
    # Else, this overwrite policy is unrecognized. Raise an exception.
    else:
        raise BetseDirException(
            'Overwrite policy "{}" unrecognized.'.format(overwrite_policy))

# ....................{ MAKERS                             }....................
@type_check
def make_unless_dir(*dirnames: str) -> None:
    '''
    Create all passed directories that do *not* already exist, silently
    ignoring those that *do* already exist.

    All nonexistent parents of this directory are also recursively created,
    reproducing the action of the POSIX-compliant ``mkdir -p`` shell command.

    Parameters
    -----------
    pathnames : tuple[str]
        Tuple of the absolute or relative pathnames of all directories to
        create.

    See Also
    -----------
    :func:`make_parent_unless_dir`
        Related function creating the parent directory of this path.
    '''

    # If this directory does *NOT* already exist, create this directory. To
    # support logging, this condition is explicitly tested for. To avoid race
    # conditions (e.g., in the event this directory is created between testing
    # and creating this directory), we preserve the makedirs() keyword argument
    # "exist_ok=True" below.
    for dirname in dirnames:
        if not is_dir(dirname):
            # Log this creation.
            logs.log_debug('Creating directory: %s', dirname)

            # Create this directory if still needed.
            os.makedirs(dirname, exist_ok=True)


@type_check
def make_parent_unless_dir(*pathnames: str) -> None:
    '''
    Create the parent directories of all passed directories that do *not*
    already exist, silently ignoring those that *do* already exist.

    Parameters
    -----------
    pathnames : tuple[str]
        Tuple of the absolute or relative pathnames to create the parent
        directories of.

    See Also
    -----------
    :func:`make_unless_dir`
        Further details.
    '''

    # Avoid circular import dependencies.
    from betse.util.path.pathnames import get_dirname, canonicalize

    # Canonicalize each pathname *BEFORE* attempting to get its dirname.
    # Relative pathnames do *NOT* have sane dirnames (e.g., the dirname for a
    # relative pathname "metatron" is the empty string) and hence *MUST* be
    # converted to absolute pathnames first.
    for pathname in pathnames:
        make_unless_dir(get_dirname(canonicalize(pathname)))

# ....................{ MAKERS ~ convenience               }....................
def canonicalize_and_make_unless_dir(dirname: str) -> str:
    '''
    Create the directory with the passed absolute or relative path if this
    directory does not already exist and return the **canonical form** (i.e.,
    unique absolute path) of this path.

    This convenience function chains the lower-level
    :func:`betse.util.path.pathnames.canonicalize` and
    :func:`make_unless_dir` functions.
    '''

    # Avoid circular import dependencies.
    from betse.util.path import pathnames

    # Dirname canonicalized from this possibly non-canonical dirname.
    dirname = pathnames.canonicalize(dirname)

    # Create this directory if needed.
    make_unless_dir(dirname)

    # Return this dirname.
    return dirname


def join_and_make_unless_dir(*partnames: str) -> str:
    '''
    Join (i.e., concatenate) the passed pathnames with the directory separator
    specific to the current platform, create the directory with the resulting
    absolute or relative path, and return this path.

    This convenience function chains the lower-level
    :func:`betse.util.path.pathnames.join` and :func:`make_unless_dir`
    functions.
    '''

    # Avoid circular import dependencies.
    from betse.util.path import pathnames

    # Dirname joined from these pathnames.
    dirname = pathnames.join(*partnames)

    # Create this directory if needed.
    make_unless_dir(dirname)

    # Return this dirname.
    return dirname

# ....................{ ITERATORS                          }....................
@type_check
def iter_basenames(dirname: str) -> SequenceTypes:
    '''
    Sequence of the basenames of all paths (e.g., non-directory files,
    subdirectories) directly contained within the passed directory.

    Parameters
    -----------
    dirname : str
        Absolute or relative dirname of the parent directory to inspect.

    Returns
    -----------
    SequenceTypes
        Sequence of the basenames of all child paths of this parent directory.
    '''

    # If this directory does *NOT* exist, raise an exception.
    die_unless_dir(dirname)

    # Sequence of all such basenames.
    return os.listdir(dirname)

# ....................{ ITERATORS ~ subdir                 }....................
@type_check
def iter_subdirnames(dirname: str) -> GeneratorType:
    '''
    Generator non-recursively yielding the dirname of each direct subdirectory
    of the passed parent directory.

    Parameters
    -----------
    dirname : str
        Absolute or relative dirname of the parent directory to inspect.

    Returns
    -----------
    GeneratorType
        Generator yielding direct subdirectory pathnames.

    Yields
    -----------
    str
        Absolute or relative dirname of each direct subdirectory of this
        directory, depending on whether the passed dirname is itself either
        absolute or relative.

    See Also
    -----------
    https://stackoverflow.com/a/25705093/2809027
        StackOverflow answer strongly inspiring this implementation.
    '''

    # Sequence of the basenames of all subdirectories of this directory
    # implemented as follows:
    #
    # * The generator returned by os.walk() is iterated once and then promptly
    #   discarded, yielding:
    #   * The ignorable absolute or relative path of this parent directory,
    #     which we (of course) were passed as input and thus already have.
    #   * This desired sequence.
    #   * An ignorable sequence of the basenames of all files of this
    #     directory.
    # * Errors emitted by low-level functions called by os.walk() (e.g.,
    #   os.listdir()) are *NOT* silently ignored.
    _, subdir_basenames, _ = next(_walk(dirname))

    # Return a generator comprehension prefixing each such basename by the
    # absolute or relative path of the parent directory.
    return (
        os_path.join(dirname, subdir_basename)
        for subdir_basename in subdir_basenames
    )


@type_check
def iter_subdir_basenames(dirname: str) -> SequenceTypes:
    '''
    Sequence of the basenames of all direct subdirectories of the passed parent
    directory.

    Parameters
    -----------
    dirname : str
        Absolute or relative dirname of the parent directory to inspect.

    Returns
    -----------
    SequenceTypes
        Sequence of direct subdirectory basenames of this parent directory.
    '''

    # Sequence of the basenames of all subdirectories of this directory. See
    # the iter_subdirnames() function for further details.
    _, subdir_basenames, _ = next(_walk(dirname))

    # Return this sequence of basenames as is.
    return subdir_basenames

# ....................{ RECURSORS ~ subdir                 }....................
@type_check
def recurse_subdirnames(dirname: str) -> GeneratorType:
    '''
    Generator recursively yielding the dirname of the passed parent directory
    and each direct *and* transitive subdirectory of that directory.

    To preserve recursive guarantees, the first dirname yielded by this
    generator is guaranteed to be the passed dirname as is.

    Parameters
    -----------
    dirname : str
        Absolute or relative dirname of the parent directory to inspect.

    Returns
    -----------
    GeneratorType
        Generator yielding direct *and* transitive subdirectory dirnames.

    Yields
    -----------
    str
        Absolute or relative dirname of each direct *and* transitive
        subdirectory of this directory, depending on whether the passed dirname
        is itself either absolute or relative.
    '''

    # Return a generator comprehension trivially yielding the absolute or
    # relative dirname of each subdirectory of this directory.
    return (subdirname for subdirname, _, _ in _walk(dirname))

# ....................{ REMOVERS                           }....................
@type_check
def remove_dir(dirname: str) -> None:
    '''
    Recursively remove the directory with the passed dirname.

    Caveats
    ----------
    Since recursive directory removal is an inherently dangerous operation,
    this function (in order):

    1. Notifies the end user with a logged warning.
    1. Waits several seconds, enabling sufficiently aware end users to jam the
       panic button.
    1. Recursively removes this directory.

    Parameters
    ----------
    dirname : str
        Absolute or relative dirname of the directory to be removed.

    Raises
    ----------
    BetseDirException
        If this directory does *not* exist.
    '''

    # Number of seconds to busywait before removing this directory.
    SLEEP_SECONDS = 4

    # Log this removal.
    logs.log_warning(
        'Removing directory in %d seconds: %s', dirname, SLEEP_SECONDS)

    # If this directory does *NOT* exist, raise an exception.
    die_unless_dir(dirname)

    # Busywait this number of seconds.
    time.sleep(SLEEP_SECONDS)

    # Recursively remove this directory.
    shutil.rmtree(dirname)

    # Log this successful completion.
    logs.log_info('Directory removed.')

# ....................{ PRIVATE ~ constants                }....................
_DIR_OVERWRITE_POLICIES_COPYTREE = frozenset((
    DirOverwritePolicy.HALT_WITH_EXCEPTION,
    DirOverwritePolicy.OVERWRITE,
))
'''
Frozen set of all **copytree-friendly directory overwrite policies** (i.e.,
:class:`.DirOverwritePolicy` enumeration members suitable for passing as the
``overwrite_policy`` parameter to the :func:`.copy_dir` function such that the
resulting implementation reduces to a trivial call to the standard
:func:`shutil.copytree` function).
'''

# ....................{ PRIVATE ~ raisers                  }....................
@type_check
def _raise_exception_dir(dirname: str) -> CallableTypes:
    '''
    Closure raising an application-specific exception encapsulating a passed
    OS exception with additional metadata describing the top-level directory
    against which this exception occurred.

    By default, exceptions raised by standard path functions (e.g.,
    :func:`os.walk`) do *not* expose this directory. When functions
    recursively change directories, these exceptions typically have
    non-human-readable messages ala:

        PermissionError: [Errno 13] Permission denied: '__pycache__'

    This closure encapsulates these low-level exceptions with higher-level
    exceptions having human-readable messages ala:

        BetsePathException: [Errno 13] Permission denied: '__pycache__' in
        "/home/leycec/py/betse/betse/util/path/".

    Design
    -----------
    This closure is principally intended to be passed as the value of the
    ``onerror`` parameter accepted by the :func:`os.walk` and :func:`os.fwalk`
    functions, preventing errors emitted by low-level functions called by these
    functions (e.g., :func:`os.listdir`) from being ignored. (By default, these
    functions silently ignore a subset of these errors.)

    Parameters
    -----------
    dirname : str
        Absolute pathname of the top-level directory to raise exceptions
        against in this closure.

    Returns
    -----------
    CallableTypes
        Closure raising such exceptions.
    '''

    @type_check
    def _raise_exception_dir_inner(exception: Exception) -> None:
        '''
        Closure encapsulating an exception passed by low-level path functions
        (e.g., :func:`os.walk`) against the absolute pathname of the top-level
        directory passed to the :func:`_raise_exception_dir` factory function.

        Parameters
        -----------
        exception : Exception
            Low-level exception raised by a low-level path function.

        Raises
        -----------
        BetsePathException
            High-level exception encapsulating this low-level exception.
        '''

        # Encapsulate this low-level exception with a higher-level exception.
        raise BetsePathException(
            '{} in "{}/".'.format(str(exception), dirname)
        ) from exception

    # Return this closure.
    return _raise_exception_dir_inner

# ....................{ PRIVATE ~ walkers                  }....................
# Undocumented os.fwalk() and os.walk() wrappers defaulting to a sane "onerror"
# callable (namely, the _raise_exception_dir() closure defined above), but
# otherwise identical to their standard variants.

def _fwalk(top, *args, **kwargs) -> GeneratorType:

    # If no "onerror" parameter was explicitly passed, default to a sane
    # exception handler.
    if 'onerror' not in kwargs:
        kwargs['onerror'] = _raise_exception_dir(top)

    # Wrap the standard variant of this generator.
    return os.fwalk(top, *args, **kwargs)


def _walk(top, *args, **kwargs) -> GeneratorType:

    # If no "onerror" parameter was explicitly passed, default to a sane
    # exception handler.
    if 'onerror' not in kwargs:
        kwargs['onerror'] = _raise_exception_dir(top)

    # Wrap the standard variant of this generator.
    return os.walk(top, *args, **kwargs)


