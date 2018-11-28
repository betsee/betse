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
from betse.util.meta.metaappabc import MetaAppABC
from betse.util.type.decorator.decmemo import property_cached
from betse.util.type.types import type_check, ModuleType

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
        return dirs.join_and_die_unless_dir(self.data_dirname, 'yaml')

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
        return files.join_and_die_unless_file(
            self.data_yaml_dirname, 'sim_config.yaml')

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

# ....................{ SINGLETONS                        }....................
#FIXME: Sadly insufficient. Downstream consumers should be trivially permitted
#to supply their own "app_meta" object overriding this default. To do so:
#
#* Preserve the assignment below. It's essential that *SOME* "app_meta"
#  attribute exist early in application startup, even if the "app_meta" object
#  in question isn't quite what the downstream consumer might expect for a
#  brief temporary period of early startup logic. (Document this, please.)
#* Define a new set_app_meta() function of this submodule with this signature:
#    def set_app_meta(app_meta: BetseMetaApp) -> None:
#  #FIXME: Actually, the above might be overkill. See the next next item.
#* Refactor the external "betsee.guimetaapp" submodule as follows:
#  * Subclass "BetseeMetaApp" from "BetseMetaApp" instead.
#  * Refactor the following line:
#    # ...from this:
#    app_meta = BetseeMetaApp()
#
#    #FIXME: Actually, this might be overkill. See the next item.
#    # ...into this:
#    from betse import metaapp
#
#    def init() -> None:
#        app_meta = BetseeMetaApp()
#        metaapp.set_app_meta(app_meta)
#* Call the betsee.guimetaapp.init() function *AS EARLY AS FEASIBLE* --
#  ideally, prior to the first call to the betse.ignition.init() function.
#  *ALTERNATELY,* we might avoid all of the above boilerplate as follows:
#  * Privatize the "betse.metaapp.app_meta" global to "_app_meta".
#  * Define a new betse.metaapp.get_app_meta() function that:
#    * Raises an exception if "_app_meta" is None.
#    * Otherwise returns "_app_meta".
#  * Refactor all references to "betse.metaapp.app_meta" to instead call
#    betse.metaapp.get_app_meta().
#  * Define a new betse.metaapp.init() function with the signature:
#      def init(app_meta: BetseMetaApp) -> None:
#          global _app_meta
#          _app_meta = app_meta
#  * Default "_app_meta" to None rather than "BetseMetaApp".
#  * Refactor the betse.ignition.init() function to have this signature:
#    def init(app_meta: BetseMetaAppOrNoneTypes) -> None:
#  * At the very start of that function:
#    * If "app_meta" is None, default that local to a new "BetseMetaApp()"
#      instance.
#    * Perform the following:
#        from betse import metaapp
#        metaapp.init(app_meta)
#* Refactor the call to the betse.ignition.init() function from the BETSEE
#  codebase to pass an instance of "BetseeMetaApp".
#* Replace all references to "betsee.guimetaapp.app_meta" with calls to
#  betse.metaapp.get_app_meta() instead; then excise the former global.
app_meta = BetseMetaApp()
'''
**Application metadata singleton** (i.e., application-wide object synopsizing
application metadata via read-only properties).
'''
