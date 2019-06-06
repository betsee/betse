#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
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
from betse import metadata, metadeps
from betse.util.app import apppath
from betse.util.app.meta.appmetaabc import AppMetaABC
from betse.util.type.decorator.decmemo import property_cached
from betse.util.type.types import type_check, ModuleType

# ....................{ SUBCLASSES                        }....................
class BetseAppMeta(AppMetaABC):
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

    # ..................{ SUPERCLASS ~ properties           }..................
    @property
    def _module_metadata(self) -> ModuleType:
        return metadata

    @property
    def _module_metadeps(self) -> ModuleType:
        return metadeps

    # ..................{ PROPERTIES ~ dir                  }..................
    @property_cached
    def betse_data_dirname(self) -> str:
        '''
        Absolute dirname of BETSE's top-level data directory if found *or*
        raise an exception otherwise (i.e., if this directory is *not* found).

        This directory typically contains BETSE-internal resources (e.g., media
        files) required at BETSE runtime.

        Design
        ----------
        The dirname returned by this property is guaranteed to be identical to
        the dirname returned by the :meth:`data_dirname` *only* if the current
        application is BETSE. If the current application is instead a
        downstream consumer of BETSE (e.g., BETSEE, PLIMBO), the dirname
        returned by this property is guaranteed to be BETSE-specific while that
        returned by the :meth:`data_dirname` property is guaranteed to be
        downstream-specific and hence differ from the former.

        In short, this property enables downstream consumers to access
        BETSE-specific data in a non-ambiguous manner *without* inviting
        inheritance or subclass issues.

        Raises
        ----------
        BetseDirException
            If this directory does *not* exist.
        '''

        # Return the absolute dirname of this BETSE-relative directory if
        # this directory exists *OR* raise an exception otherwise.
        return apppath.get_dirname(package=betse, dirname='data')


    @property_cached
    def betse_data_yaml_dirname(self) -> str:
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
        return dirs.join_or_die(self.betse_data_dirname, 'yaml')

    # ..................{ PROPERTIES ~ file                 }..................
    @property_cached
    def betse_sim_conf_default_filename(self) -> str:
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
        return files.join_or_die(
            self.betse_data_yaml_dirname, 'sim_config.yaml')

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
