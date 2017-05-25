#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
High-level constants describing this application's filesystem usage.

These constants provide the absolute paths of files and directories intended for
general use by both this application and downstream reverse dependencies of this
application (e.g., BETSEE, the BETSE GUI). For portability, these constants are
initialized in a system-aware manner guaranteed to be sane under various
installation environments -- including PyInstaller-frozen executables and
:mod:`setuptools`-installed script wrappers.
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
from betse.util.type.call.memoizers import callable_cached
from betse.util.type.types import type_check

# ....................{ GETTERS ~ dir                      }....................
@callable_cached
def get_root_dirname() -> str:
    '''
    Absolute path of the root directory suffixed by a directory separator.

    The definition of "root directory" conditionally depends on the current
    platform. If this platform is:

    * POSIX-compatible (e.g., Linux, OS X), this is simply ``/``.
    * Microsoft Windows, this is the value of the ``%HOMEDRIVE%`` environment
      variable. This is the ``:``-suffixed letter of the drive to which Windows
      was originally installed -- typically but *not* necessarily ``C:\``.
    '''

    # Avoid circular import dependencies.
    from betse.util.os import oses
    from betse.util.os.shell import shellenv

    # Return this dirname.
    if oses.is_windows_vanilla():
        return shellenv.get_var_or_default('HOMEDRIVE', 'C:') + os.path.sep
    else:
        return os.path.sep


@callable_cached
def get_home_dirname() -> str:
    '''
    Absolute path of the home directory of the current user if found *or* raise
    an exception otherwise (i.e., if this directory is *not* found).
    '''

    # Avoid circular import dependencies.
    from betse.util.path import dirs, pathnames

    # Absolute path of this directory.
    home_dirname = pathnames.canonicalize('~')

    # If this directory is not found, fail.
    dirs.die_unless_dir(home_dirname)

    # Return this directory's path.
    return home_dirname


@callable_cached
def get_dot_dirname() -> str:
    '''
    Absolute path of this application's top-level dot directory in the home
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
@callable_cached
def get_data_dirname() -> str:
    '''
    Absolute path of this application's top-level data directory if found *or*
    raise an exception otherwise (i.e., if this directory is *not* found).

    This directory contains application-internal resources (e.g., media files)
    required at application runtime.
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


@callable_cached
def get_data_yaml_dirname() -> str:
    '''
    Absolute path of this application's data subdirectory containing
    YAML-formatted files if found *or* raise an exception otherwise (i.e., if
    this directory is *not* found).

    This directory contains:

    * The YAML-formatted plaintext file specifying a simulation configuration.
    * All assets referenced in and hence required by this file.
    '''

    # Avoid circular import dependencies.
    from betse.util.path import dirs

    # Return this dirname if this directory exists or raise an exception.
    return dirs.join_and_die_unless_dir(get_data_dirname(), 'yaml')


@callable_cached
def get_data_default_asset_dirnames() -> str:
    '''
    Tuple of the absolute paths of all directories containing assets (e.g.,
    media files, YAML files) referenced and hence required by this application's
    default simulation configuration file if found *or* raise an exception
    otherwise (i.e., if any such directory is *not* found).

    This includes directories containing (in no particular order):

    * YAML-formatted configuration files specifying fine-grained simulation
      settings, including:
      * The biochemical reaction network (BRN) to be simulated if any.
      * The gene regulatory network (GRN) to be simulated if any.
    * Ancillary media files, including:
      * Sample images typically used as input geometries for tissue profiles.
        Each image defines a cell cluster shape typically corresponding to some
        species (e.g., planaria) or organ (e.g., heart), complete with internal
        cellular structure and exterior cellular boundary. Each spatial
        subdivision of this shape is then associated with a real-world tissue
        profile in the configuration file(s) referencing this image.
    '''

    # Avoid circular import dependencies.
    from betse.util.path import dirs, pathnames

    # Absolute path of the top-level YAML-specific data directory.
    DATA_YAML_DIRNAME = get_data_yaml_dirname()

    # Tuple of the absolute paths of all default asset directories.
    data_default_asset_dirnames = (
        pathnames.join(DATA_YAML_DIRNAME, 'geo'),
        pathnames.join(DATA_YAML_DIRNAME, 'extra_configs'),
    )

    # If any such directory is not found, fail.
    dirs.die_unless_dir(*data_default_asset_dirnames)

    # Return this tuple.
    return data_default_asset_dirnames

# ....................{ GETTERS ~ file                     }....................
@callable_cached
def get_log_default_filename() -> str:
    '''
    Absolute path of this application's default user-specific logfile.

    This is the plaintext file to which all messages are logged by default.
    '''

    # Avoid circular import dependencies.
    from betse.util.path import pathnames

    # Return the absolute path of this file.
    return pathnames.join(get_dot_dirname(), metadata.SCRIPT_BASENAME + '.log')


@callable_cached
def get_profile_default_filename() -> str:
    '''
    Absolute path of this application's default user-specific profile dumpfile.

    This is the binary file to which profiled statistics are saved by default.
    '''

    # Avoid circular import dependencies.
    from betse.util.path import pathnames

    # Return the absolute path of this file.
    return pathnames.join(get_dot_dirname(), metadata.SCRIPT_BASENAME + '.prof')


@callable_cached
def get_sim_config_default_filename() -> str:
    '''
    Absolute path of this application's default simulation configuration file
    if found *or* raise an exception otherwise (i.e., if this file is *not*
    found).

    This is the plaintext YAML-formatted file from which all new simulation
    configurations derive.
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
    Absolute path of this application's default user-specific history file for
    the read–eval–print loop (REPL) implemented by the third-party module with
    the passed name.

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
