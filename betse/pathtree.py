#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Collection of the absolute paths of numerous critical files and directories
describing the structure of this application on the local filesystem.

These are intended for consumption by both this application and downstream
reverse dependencies of this application (e.g., BETSEE, the BETSE GUI). For
portability, these paths are initialized in a system-aware manner guaranteed to
be sane under insane installation environments -- including PyInstaller-frozen
executables and :mod:`setuptools`-installed script wrappers.
'''

#FIXME: The current globals-based approach is inefficient in the case of BETSE
#being installed as a compressed EGG rather than an uncompressed directory. In
#the former case, the current approach (namely, the call to
#resources.get_pathname() performed below) silently extracts the entirety of
#this egg to a temporary setuptools-specific cache directory. That's bad. To
#circumvent this, we'll need to refactor the codebase to directly require only
#"file"-like objects rather than indirectly requiring the absolute paths of
#data resources that are then opened as "file"-like objects.
#
#Specifically, whenever we require a "file"-like object for a codebase resource,
#we'll need to call the setuptools-specific pkg_resources.resource_stream()
#function rather than attempting to open the path given by a global below.
#Ultimately, *ALL* of the codebase-specific globals declared below (e.g.,
#"DATA_DIRNAME") should go away.

# ....................{ IMPORTS                            }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To raise human-readable exceptions on missing mandatory dependencies,
# the top-level of this module may import *ONLY* from packages guaranteed to
# exist at installation time (i.e., standard Python packages). Likewise, to
# avoid circular import dependencies, the top-level of this module should avoid
# importing application packages except where explicitly required.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
import os
from betse import metadata
from betse.exceptions import BetseModuleException
from betse.util.type.decorator.decmemo import func_cached
from betse.util.type.types import type_check, StrOrNoneTypes

# ....................{ GETTERS ~ dir                      }....................
@func_cached
def get_root_dirname() -> str:
    '''
    Absolute pathname of the root directory suffixed by a directory separator.

    The definition of "root directory" conditionally depends on the current
    platform. If this platform is:

    * POSIX-compatible (e.g., Linux, OS X), this is simply ``/``.
    * Microsoft Windows, this is the value of the ``%HOMEDRIVE%`` environment
      variable. This is the ``:``-suffixed letter of the drive to which Windows
      was originally installed -- typically but *not* necessarily ``C:\\``.
    '''

    # Avoid circular import dependencies.
    from betse.util.os import oses
    from betse.util.os.shell import shellenv

    # Return this dirname.
    if oses.is_windows_vanilla():
        return shellenv.get_var_or_default('HOMEDRIVE', 'C:') + os.path.sep
    else:
        return os.path.sep


@func_cached
def get_home_dirname() -> str:
    '''
    Absolute pathname of the home directory of the current user if found *or*
    raise an exception otherwise (i.e., if this directory is *not* found).
    '''

    # Avoid circular import dependencies.
    from betse.util.path import dirs, pathnames

    # Absolute path of this directory.
    home_dirname = pathnames.canonicalize('~')

    # If this directory is not found, fail.
    dirs.die_unless_dir(home_dirname)

    # Return this directory's path.
    return home_dirname

# ....................{ GETTERS ~ dir : app                }....................
@func_cached
def get_dot_dirname() -> str:
    '''
    Absolute pathname of this application's top-level dot directory in the home
    directory of the current user, silently creating this directory if *not*
    already found.

    This directory contains user-specific files (e.g., logfiles, profile files)
    both read from and written to at application runtime. These are typically
    plaintext files consumable by external users and third-party utilities.

    Locations
    ----------
    Specifically, this path is:

    * Under Linux, ``~/.betse/``. BETSE does *not* currently comply with the
      _XDG Base Directory Specification (e.g., ``~/.local/share/betse``), which
      the principal authors of BETSE regard as unuseful (if not harmful).
    * Under OS X, ``~/Library/Application Support/betse``.
    * Under Windows,
      ``C:\\Documents and Settings\\${User}\\Application Data\\betse``.

    .. _XDG Base Directory Specification:
        http://standards.freedesktop.org/basedir-spec/basedir-spec-latest.html

    Raises
    ----------
    BetseDirException
        If this directory does *not* exist.
    '''

    # Avoid circular import dependencies.
    from betse.util.os import oses
    from betse.util.os.shell import shellenv
    from betse.util.path import dirs, pathnames

    # Absolute path of this directory.
    dot_dirname = None

    # If the current platform is macOS, return the appropriate directory.
    if oses.is_macos():
        dot_dirname = pathnames.join(
            get_home_dirname(),
            'Library',
            'Application Support',
            metadata.SCRIPT_BASENAME,
        )
    # If the current platform is Windows, return the appropriate directory.
    elif oses.is_windows():
        dot_dirname = pathnames.join(shellenv.get_var('APPDATA'), metadata.NAME)
    # Else, assume the current platform to be POSIX-compatible.
    else:
        #FIXME: Explicitly assert POSIX compatibility here. To do so, we'll want
        #to define and call a new betse.util.os.oses.die_unless_posix()
        #function here.

        dot_dirname = pathnames.join(
            get_home_dirname(), '.' + metadata.SCRIPT_BASENAME)

    # Create this directory if not found.
    dirs.make_unless_dir(dot_dirname)

    # Return this directory's path.
    return dot_dirname

# ....................{ GETTERS ~ dir : data               }....................
@func_cached
def get_data_dirname() -> str:
    '''
    Absolute pathname of this application's top-level data directory if found
    *or* raise an exception otherwise (i.e., if this directory is *not* found).

    This directory contains application-internal resources (e.g., media files)
    required at application runtime.

    Raises
    ----------
    BetseDirException
        If this directory does *not* exist.
    '''

    # Avoid circular import dependencies.
    import betse
    from betse.util.path import dirs, pathnames

    # Absolute path of this directory.
    data_dirname = pathnames.get_app_pathname(package=betse, pathname='data')

    # If this directory is not found, raise an exception.
    dirs.die_unless_dir(data_dirname)

    # Return the absolute path of this directory.
    return data_dirname


@func_cached
def get_data_yaml_dirname() -> str:
    '''
    Absolute pathname of this application's data subdirectory containing
    YAML-formatted files if found *or* raise an exception otherwise (i.e., if
    this directory is *not* found).

    This directory contains:

    * The YAML-formatted plaintext file specifying a simulation configuration.
    * All assets referenced in and hence required by this file.

    Raises
    ----------
    BetseDirException
        If this directory does *not* exist.
    '''

    # Avoid circular import dependencies.
    from betse.util.path import dirs

    # Return this dirname if this directory exists or raise an exception.
    return dirs.join_and_die_unless_dir(get_data_dirname(), 'yaml')

# ....................{ GETTERS ~ dir : package            }....................
@func_cached
def get_package_dirname() -> str:
    '''
    Absolute pathname of this application's top-level package directory if found
    *or* raise an exception otherwise (i.e., if this directory is *not* found).

    This directory typically resides in the ``site-packages`` subdirectory of
    the system-wide standard lib for the active Python interpreter (e.g.,
    ``/usr/lib64/python3.6/site-packages/betse``).

    Raises
    ----------
    BetseDirException
        If this directory does *not* exist.
    '''

    # Avoid circular import dependencies.
    import betse
    from betse.util.path import dirs
    from betse.util.type import modules

    # Absolute pathname of the directory yielding the top-level "betse" package.
    package_dirname = modules.get_dirname(betse)

    # If this directory is not found, fail.
    dirs.die_unless_dir(package_dirname)

    # Return this directory's pathname.
    return package_dirname

# ....................{ GETTERS ~ dir : git                }....................
@func_cached
def get_git_worktree_dirname() -> str:
    '''
    Absolute pathname of this application's Git-based **working tree** (i.e.,
    top-level directory containing this application's ``.git`` subdirectory and
    ``setup.py`` install script) if this application was installed in a
    developer manner *or* raise an exception otherwise (i.e., if this directory
    is *not* found).

    Raises
    ----------
    BetseDirException
        If this directory does *not* exist.
    '''

    # Avoid circular import dependencies.
    from betse.util.path import dirs

    # Absolute pathname of this application's Git-based working tree if this
    # application was installed in a developer manner or "None" otherwise.
    git_worktree_dirname = get_git_worktree_dirname_or_none()

    # If this directory is not found, fail.
    dirs.die_unless_dir(git_worktree_dirname)

    # Return this directory's pathname.
    return git_worktree_dirname


@func_cached
def get_git_worktree_dirname_or_none() -> StrOrNoneTypes:
    '''
    Absolute pathname of this application's Git-based **working tree** (i.e.,
    top-level directory containing this application's ``.git`` subdirectory and
    ``setup.py`` install script) if this application was installed in a
    developer manner *or* ``None`` otherwise.

    Returns
    ----------
    StrOrNoneTypes
        Specifically, this function returns either:
        * A pathname if this application was installed in a developer manner,
          typically either via:
          * ``python3 setup.py develop``.
          * ``python3 setup.py symlink``.
        * ``None`` if this application was installed in a non-developer manner,
          typically either via:
          * ``pip3 install``.
          * ``python3 setup.py develop``.
    '''

    # Avoid circular import dependencies.
    from betse.util.path import dirs, pathnames

    # Absolute pathname of the directory offering the top-level "betse" package,
    # canonicalized into a directory rather than symbolic link to increase the
    # likelihood of obtaining the actual parent directory of this package.
    package_dirname = pathnames.canonicalize(get_package_dirname())

    # Absolute pathname of the parent directory of this directory.
    worktree_dirname = pathnames.get_dirname(package_dirname)

    # Absolute pathname of the ".git" subdirectory of this parent directory.
    git_subdirname = pathnames.join(worktree_dirname, '.git')

    # Return this parent directory's absolute pathname if this subdirectory
    # exists *OR* "None" otherwise.
    return worktree_dirname if dirs.is_dir(git_subdirname) else None

# ....................{ GETTERS ~ file                     }....................
@func_cached
def get_log_default_filename() -> str:
    '''
    Absolute pathname of this application's default user-specific logfile.

    This is the plaintext file to which all messages are logged by default.
    '''

    # Avoid circular import dependencies.
    from betse.util.path import pathnames

    # Return the absolute path of this file.
    return pathnames.join(get_dot_dirname(), metadata.SCRIPT_BASENAME + '.log')


@func_cached
def get_profile_default_filename() -> str:
    '''
    Absolute pathname of this application's default user-specific profile
    dumpfile.

    This is the binary file to which profiled statistics are saved by default.
    '''

    # Avoid circular import dependencies.
    from betse.util.path import pathnames

    # Return the absolute path of this file.
    return pathnames.join(get_dot_dirname(), metadata.SCRIPT_BASENAME + '.prof')


@func_cached
def get_sim_config_default_filename() -> str:
    '''
    Absolute pathname of this application's default simulation configuration
    file if found *or* raise an exception otherwise (i.e., if this file is *not*
    found).

    This is the plaintext YAML-formatted file from which all new simulation
    configurations derive.

    Raises
    ----------
    BetseFileException
        If this file does *not* exist.
    '''

    # Avoid circular import dependencies.
    from betse.util.path import files

    # Return this dirname if this directory exists or raise an exception.
    return files.join_and_die_unless_file(
        get_data_yaml_dirname(), 'sim_config.yaml')

# ....................{ GETTERS ~ file : arg               }....................
@type_check
def get_repl_history_filename(repl_module_name: str) -> dict:
    '''
    Absolute pathname of this application's default user-specific history file
    for the read–eval–print loop (REPL) implemented by the third-party module
    with the passed name.

    This history file is a plaintext file to which this REPL appends each read
    user command for preserving command history between REPL sessions.

    Parameters
    ----------
    repl_module_name : str
        Fully-qualified name of the third-party module implementing this REPL.

    Returns
    ----------
    str
        Absolute path of the user-specific history file for this REPL.

    Raises
    ----------
    BetseModuleException
        If this REPL module name is unrecognized.
    '''

    # Avoid circular import dependencies.
    from betse.util.path import pathnames

    # Absolute path of this application's user-specific dot directory.
    DOT_DIRNAME = get_dot_dirname()

    # Dictionary mapping each REPL module name to its history filename.
    REPL_MODULE_NAME_TO_HISTORY_FILENAME = {
        'ptpython': pathnames.join(DOT_DIRNAME, 'ptpython.hist'),
        'readline': pathnames.join(DOT_DIRNAME, 'readline.hist'),
    }

    # If this REPL module name is unrecognized, fail.
    if repl_module_name not in REPL_MODULE_NAME_TO_HISTORY_FILENAME:
        raise BetseModuleException(
            'REPL module "{}" unrecognized.'.format(repl_module_name))

    # Else, return the history filename for this REPL.
    return REPL_MODULE_NAME_TO_HISTORY_FILENAME[repl_module_name]
