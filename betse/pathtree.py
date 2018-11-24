#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
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

#FIXME: Refactor the remainder of this submodule into the newly annointed
#"betse.util.meta.metaappabc" submodule.

# ....................{ IMPORTS                           }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To raise human-readable exceptions on missing mandatory dependencies,
# the top-level of this module may import *ONLY* from packages guaranteed to
# exist at installation time (i.e., standard Python packages). Likewise, to
# avoid circular import dependencies, the top-level of this module should avoid
# importing application packages except where explicitly required.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

from betse import metadata
from betse.metaapp import app_meta
from betse.exceptions import BetseModuleException
from betse.util.type.decorator.decmemo import func_cached
from betse.util.type.types import type_check

# ....................{ GETTERS ~ dir : data              }....................
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
    return dirs.join_and_die_unless_dir(app_meta.data_dirname, 'yaml')

# ....................{ GETTERS ~ dir : package           }....................
@func_cached
def get_package_dirname() -> str:
    '''
    Absolute pathname of this application's top-level package directory if
    found *or* raise an exception otherwise (i.e., if this directory is *not*
    found).

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
    from betse.util.py import pymodule

    # Absolute pathname of the directory yielding the top-level "betse" package.
    package_dirname = pymodule.get_dirname(betse)

    # If this directory is not found, fail.
    dirs.die_unless_dir(package_dirname)

    # Return this directory's pathname.
    return package_dirname

# ....................{ GETTERS ~ file                    }....................
@func_cached
def get_log_default_filename() -> str:
    '''
    Absolute filename of this application's default user-specific logfile.

    This is the plaintext file to which all messages are logged by default.
    '''

    # Avoid circular import dependencies.
    from betse.util.path import pathnames

    # Return the absolute path of this file.
    return pathnames.join(
        app_meta.dot_dirname, metadata.SCRIPT_BASENAME + '.log')


@func_cached
def get_profile_default_filename() -> str:
    '''
    Absolute filename of this application's default user-specific profile
    dumpfile.

    This is the binary file to which profiled statistics are saved by default.
    '''

    # Avoid circular import dependencies.
    from betse.util.path import pathnames

    # Return the absolute path of this file.
    return pathnames.join(
        app_meta.dot_dirname, metadata.SCRIPT_BASENAME + '.prof')


@func_cached
def get_sim_config_default_filename() -> str:
    '''
    Absolute pathname of this application's default simulation configuration
    file if found *or* raise an exception otherwise (i.e., if this file is
    *not* found).

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

# ....................{ GETTERS ~ file : arg              }....................
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

    # Dictionary mapping each REPL module name to its history filename.
    REPL_MODULE_NAME_TO_HISTORY_FILENAME = {
        'ptpython': pathnames.join(app_meta.dot_dirname, 'ptpython.hist'),
        'readline': pathnames.join(app_meta.dot_dirname, 'readline.hist'),
    }

    # If this REPL module name is unrecognized, fail.
    if repl_module_name not in REPL_MODULE_NAME_TO_HISTORY_FILENAME:
        raise BetseModuleException(
            'REPL module "{}" unrecognized.'.format(repl_module_name))

    # Else, return the history filename for this REPL.
    return REPL_MODULE_NAME_TO_HISTORY_FILENAME[repl_module_name]
