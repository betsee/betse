#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
High-level constants describing this application's filesystem usage.

These constants provide the absolute paths of files and directories intended for
general use by both the CLI and GUI. For portability, such constants are
initialized in a system-aware manner guaranteed to be sane under various
installation environments -- including PyInstaller-frozen executables and
`setuptools`-installed script wrappers.
'''

#FIXME: Refactor all globals defined below into @callable_cached-decorated
#module functions. Note that we'll need to preserve the "HOME_DIRECTORY" global
#for backward compatibility purposes. All other globals *SHOULD* be safely
#refactorable and then removable.

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
# exist at installation time (i.e., stock Python and BETSE packages).
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
from betse import metadata
from betse.exceptions import BetseModuleException
from betse.util.type.call.memoizers import callable_cached
from betse.util.type.types import type_check
from os import environ, path

# ....................{ GETTERS ~ dir                      }....................
@callable_cached
def get_home_dirname() -> str:
    '''
    Absolute path of the home directory of the current user, raising an
    exception if this directory is *not* found.
    '''

    # Avoid circular import dependencies.
    from betse.util.path import dirs

    # Absolute path of this directory.
    home_dirname = path.expanduser('~')

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

    This directory contains application-external resources (e.g., configuration
    files) created at application runtime and subsequently editable by external
    users and utilities.

    Locations
    ----------
    Specifically, this path is:

    * Under Linux, ``~/.betse/``. BETSE does *not* currently comply with the
      _XDG Base Directory Specification (e.g., ``~/.local/share/betse``), which
      the principal authors of BETSE regard as unhelpful -- if not harmful.
    * Under OS X, ``~/Library/Application Support/betse``.
    * Under Windows,
      ``C:\Documents and Settings\${User}\Application Data\betse``.

    .. _XDG Base Directory Specification: http://standards.freedesktop.org/basedir-spec/basedir-spec-latest.html
    '''

    # Avoid circular import dependencies.
    from betse.util.path import dirs, paths
    from betse.util.os import oses

    # Absolute path of this directory.
    dot_dirname = None

    # If the current platform is macOS, return the appropriate directory.
    if oses.is_macos():
        dot_dirname = paths.join(
            get_home_dirname(),
            'Library',
            'Application Support',
            metadata.SCRIPT_NAME_CLI,
        )
    # If the current platform is Windows, return the appropriate directory.
    elif oses.is_windows():
        dot_dirname = paths.join(environ['APPDATA'], metadata.NAME)
    # Else, assume the current platform to be POSIX-compatible.
    else:
        #FIXME: Explicitly assert POSIX compatibility here. To do so, we'll want
        #to define and call a new betse.util.os.oses.die_unless_posix()
        #function here.

        dot_dirname = paths.join(
            get_home_dirname(), '.' + metadata.SCRIPT_NAME_CLI)

    # Create this directory if not found.
    dirs.make_unless_dir(dot_dirname)

    # Return this directory's path.
    return dot_dirname

# ....................{ GETTERS ~ dir : data               }....................
@callable_cached
def get_data_dirname() -> str:
    '''
    Absolute path of this application's top-level data directory, raising an
    exception if this directory is *not* found.

    This directory contains application-internal resources (e.g., media files)
    required at application runtime.
    '''

    # Avoid circular import dependencies.
    from betse.lib.setuptools import resources
    from betse.util.path import dirs, paths
    from betse.util.py import freezers

    # Absolute path of this directory.
    data_dirname = None

    # Basename of this directory relative to the directory containing this
    # submodule. Since setuptools-specific resource pathnames expect the POSIX-
    # rather than Windows-specific directory separator (i.e., "/" rather than
    # "\"), this basename must *NOT* contain the latter.
    data_basename = 'data'

    # If the current application is a PyInstaller-frozen executable binary,
    # defer to the PyInstaller-specific private attribute "_MEIPASS" added
    # to the canonical "sys" module by the PyInstaller bootloader embedded in
    # this binary. This attribute provides the absolute path of the temporary
    # directory containing all application data resources extracted from this
    # binary by this bootloader. "And it's turtles all the way down."
    if freezers.is_frozen_pyinstaller():
        data_dirname = paths.join(
            freezers.get_app_dirname_pyinstaller(), data_basename)
    # If the current application is a setuptools-installed script wrapper, the
    # data directory will have been preserved as is in the setuptools-installed
    # copy of the current Python package tree. In this case, query setuptools to
    # obtain this directory's path in a cross-platform manner.
    elif resources.is_dir(__name__, data_basename):
        data_dirname = resources.get_pathname(__name__, data_basename)
    # Else, the current application is either a setuptools-symlinked script
    # wrapper *OR* was invoked via the hidden "python3 -m betse.cli.cli"
    # command. In either case, such directory's path is directly obtainable
    # relative to the absolute path of the current module.
    else:
        data_dirname = paths.join(paths.get_dirname(__file__), data_basename)

    # If this directory is not found, fail.
    dirs.die_unless_dir(data_dirname)

    # Return this directory's path.
    return data_dirname


@callable_cached
def get_data_yaml_dirname() -> str:
    '''
    Absolute path of this application's top-level YAML-specific data directory,
    raising an exception if this directory is *not* found.

    This directory contains:

    * The YAML-formatted plaintext file specifying a simulation configuration.
    * All assets referenced in and hence required by this file.
    '''

    # Avoid circular import dependencies.
    from betse.util.path import dirs, paths

    # Absolute path of this directory.
    data_yaml_dirname = paths.join(get_data_dirname(), 'yaml')

    # If this directory is not found, fail.
    dirs.die_unless_dir(data_yaml_dirname)

    # Return this directory's path.
    return data_yaml_dirname


@callable_cached
def get_data_default_asset_dirnames() -> str:
    '''
    Tuple of the absolute paths of all directories containing assets (e.g.,
    media files, YAML files) referenced and hence required by this application's
    default simulation configuration file, raising an exception if any such
    directory is *not* found.

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
    from betse.util.path import dirs, paths

    # Absolute path of the top-level YAML-specific data directory.
    DATA_YAML_DIRNAME = get_data_yaml_dirname()

    # Tuple of the absolute paths of all default asset directories.
    data_default_asset_dirnames = (
        paths.join(DATA_YAML_DIRNAME, 'geo'),
        paths.join(DATA_YAML_DIRNAME, 'extra_configs'),
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
    from betse.util.path import paths

    # Return the absolute path of this file.
    return paths.join(get_dot_dirname(), metadata.SCRIPT_NAME_CLI + '.log')


@callable_cached
def get_profile_default_filename() -> str:
    '''
    Absolute path of this application's default user-specific profile dumpfile.

    This is the binary file to which profiled statistics are saved by default.
    '''

    # Avoid circular import dependencies.
    from betse.util.path import paths

    # Return the absolute path of this file.
    return paths.join(get_dot_dirname(), metadata.SCRIPT_NAME_CLI + '.prof')


@callable_cached
def get_sim_config_default_filename() -> str:
    '''
    Absolute path of this application's default simulation configuration file,
    raising an exception if this file does *not* exist.

    This is the plaintext YAML-formatted file from which all new simulation
    configurations derive.
    '''

    # Avoid circular import dependencies.
    from betse.util.path import files, paths

    # Absolute path of this file.
    config_default_filename = paths.join(
        get_data_yaml_dirname(), 'sim_config.yaml')

    # If this file is not found, fail.
    files.die_unless_file(config_default_filename)

    # Return this file's path.
    return config_default_filename

# ....................{ GETTERS ~ file : arg               }....................
@type_check
def get_repl_history_filename(repl_module_name: str) -> dict:
    '''
    Absolute path of this application's default user-specific history file for
    the read–eval–print loop (REPL) implemented by the third-party module with
    the passed name..

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
    from betse.util.path import paths

    # Absolute path of this application's user-specific dot directory.
    DOT_DIRNAME = get_dot_dirname()

    # Dictionary mapping each REPL module name to its history filename.
    REPL_MODULE_NAME_TO_HISTORY_FILENAME = {
        'ptpython': paths.join(DOT_DIRNAME, 'ptpython.hist'),
        'readline': paths.join(DOT_DIRNAME, 'readline.hist'),
    }

    # If this REPL module name is unrecognized, fail.
    if repl_module_name not in REPL_MODULE_NAME_TO_HISTORY_FILENAME:
        raise BetseModuleException(
            'REPL module "{}" unrecognized.'.format(repl_module_name))

    # Else, return the history filename for this REPL.
    return REPL_MODULE_NAME_TO_HISTORY_FILENAME[repl_module_name]
