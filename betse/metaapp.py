#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level **application metadata singleton** (i.e., application-wide object
synopsizing application metadata via read-only properties).
'''

# ....................{ IMPORTS                           }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: This subclass must *NOT* be accessed from the top-level "setup.py"
# script of this or any other application. This application and hence this
# subclass is *NOT* guaranteed to exist at setuptools-based installation-time
# for downstream consumers (e.g., BETSEE).
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

import betse
from betse.exceptions import BetseMetaAppException
from betse.util.meta.metaappabc import MetaAppABC
from betse.util.type.decorator.decmemo import property_cached
from betse.util.type.types import type_check, ModuleType, NoneType

# ....................{ SUBCLASSES                        }....................
class BetseMetaApp(MetaAppABC):
    '''
    **Application metadata singleton** (i.e., application-wide object
    synopsizing application metadata via read-only properties) subclass.

    Caveats
    ----------
    **This subclass must not be accessed from the top-level ``setup.py`` script
    of this or any other application.** This application and hence this
    subclass is *not* guaranteed to exist at setuptools-based installation-time
    for downstream consumers (e.g., BETSEE).
    '''

    # ..................{ PROPERTIES ~ public : superclass  }..................
    # Abstract read-only properties required to be defined by subclasses.

    @property
    def package(self) -> ModuleType:
        return betse

    # ..................{ PROPERTIES ~ dir                  }..................
    @property_cached
    def data_yaml_dirname(self) -> str:
        '''
        Absolute pathname of this application's data subdirectory containing
        YAML-formatted files if found *or* raise an exception otherwise (i.e.,
        if this directory is *not* found).

        This directory contains:

        * The YAML-formatted plaintext file specifying a simulation
          configuration.
        * All assets referenced in and hence required by this file.

        Raises
        ----------
        BetseDirException
            If this directory does *not* exist.
        '''

        # Avoid circular import dependencies.
        from betse.util.path import dirs

        # Return this dirname if this directory exists or raise an exception.
        return dirs.join_or_die(self.data_dirname, 'yaml')

    # ..................{ PROPERTIES ~ file                 }..................
    @property_cached
    def sim_conf_default_filename(self) -> str:
        '''
        Absolute pathname of this application's **default simulation
        configuration file** (i.e., plaintext YAML-formatted file from which
        all simulation configurations derive) if found *or* raise an exception
        otherwise (i.e., if this file is *not* found).

        Raises
        ----------
        BetseFileException
            If this file does *not* exist.
        '''

        # Avoid circular import dependencies.
        from betse.util.path import files

        # Return this dirname if this directory exists or raise an exception.
        return files.join_or_die(self.data_yaml_dirname, 'sim_config.yaml')

    # ..................{ GETTERS                           }..................
    @type_check
    def get_repl_history_filename(self, repl_module_name: str) -> dict:
        '''
        Absolute pathname of this application's default user-specific history
        file for the read–eval–print loop (REPL) implemented by the third-party
        module with the passed name.

        This history file is a plaintext file to which this REPL appends each
        read user command for preserving command history between REPL sessions.

        Parameters
        ----------
        repl_module_name : str
            Fully-qualified name of the third-party module implementing this
            REPL.

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
        from betse.exceptions import BetseModuleException
        from betse.util.path import pathnames

        # Dictionary mapping each REPL module name to its history filename.
        REPL_MODULE_NAME_TO_HISTORY_FILENAME = {
            'ptpython': pathnames.join(self.dot_dirname, 'ptpython.hist'),
            'readline': pathnames.join(self.dot_dirname, 'readline.hist'),
        }

        # If this REPL module name is unrecognized, fail.
        if repl_module_name not in REPL_MODULE_NAME_TO_HISTORY_FILENAME:
            raise BetseModuleException(
                'REPL module "{}" unrecognized.'.format(repl_module_name))

        # Else, return the history filename for this REPL.
        return REPL_MODULE_NAME_TO_HISTORY_FILENAME[repl_module_name]

# ....................{ TYPES                             }....................
BetseMetaAppOrNoneTypes = (BetseMetaApp, NoneType)
'''
Tuple of the types of both the application metadata *and* ``None``  singletons.
'''

# ....................{ GLOBALS                           }....................
_app_meta = None
'''
**Application metadata singleton** (i.e., application-wide object synopsizing
application metadata via read-only properties).

Caveats
----------
For safety, callers are advised to call the :func:`get_app_meta` getter safely
returning this private singleton rather than directly accessing this private
singleton unsafely.
'''

# ....................{ GETTERS                           }....................
def get_app_meta() -> BetseMetaApp:
    '''
    **Application metadata singleton** (i.e., application-wide object
    synopsizing application metadata via read-only properties) if this
    singleton has already been instantiated by a prior call to the :func:`init`
    function *or* raise an exception otherwise (i.eg., if that function has yet
    to be called).

    Returns
    ----------
    BetseMetaApp
        Application metadata singleton defined by the most recent call to the
        :func:`init` function.

    Raises
    ----------
    BetseMetaAppException
        If the :func:`init` function has yet to be called.
    '''

    # If no application metadata singleton exists, raise an exception.
    if not _app_meta:
        raise BetseMetaAppException(
            'Application metadata singleton undefined '
            '(e.g., as betse.metaapp.init() not called).')

    # Else, an application metadata singleton exists; return it, please.
    return _app_meta

# ....................{ INITIALIZERS                      }....................
@type_check
def init(app_meta: BetseMetaAppOrNoneTypes = None) -> None:
    '''
    Initialize this submodule with either the passed application metadata
    singleton if non-``None`` *or* a new instance of the :class:`BetseMetaApp`
    subclass otherwise (i.e., if no such singleton is passed).

    Parameters
    ----------
    app_meta : BetseMetaAppOrNoneTypes
        Caller-specific application metadata singleton (i.e., instance of the
        :class:`BetseMetaApp` subclass). Defaults to ``None``, in which case
        this parameter defaults to a vanilla instance of that subclass.
    '''

    # Enable this singleton global to be overwritten be the passed parameter.
    global _app_meta

    # If passed no caller-specific singleton, default to a generic singleton.
    if app_meta is None:
        app_meta = BetseMetaApp()

    # Set this singleton global to this caller-specific singleton.
    _app_meta = app_meta
