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
from betse.exceptions import BetseMetaAppException
from betse.util.app import apppath
from betse.util.app.meta.metaappabc import MetaAppABC
from betse.util.type.decorator.decmemo import property_cached
from betse.util.type.types import type_check, NoneType

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

    # ..................{ PROPERTIES ~ dir                  }..................
    #FIXME: Fundamentally broken. The issue, of course, is the access of
    #"self.data_dirname" below -- which, in subclasses, mistakenly refers to
    #the subclass data directory unlikely to contain a "yaml" subdirectory.
    #The solution is... well, it's a bit awkward. It also appears to be the
    #only reasonable means of resolving this issue. We want to:
    #
    #* Unconditionally instantiate an instance of this "BetseMetaApp" subclass
    #  in this module as a new private global singleton -- say,
    #  "_app_meta_betse = BetseMetaApp()" defined either at the top-level or
    #  (probably ideally, for safety) in the init() method.
    #* Refactor the "self.data_dirname" below to read
    #  "_app_meta_betse.data_dirname" instead.
    #
    #Naturally, subclasses may still elect to override this sane default
    #implementation if they choose -- but they really shouldn't, ever.
    #FIXME: Perhaps not? The above approach is patently absurd and blatantly
    #overkill; simply extract these properties into @callable_cached-decorated
    #public top-level functions of this submodule as under the prior design. On
    #doing so, of course, note *EXACTLY* why this is being done.
    #FIXME: Perhaps, actually. The prior commentary would be ideal, except for
    #the obvious conundrum of needing to access the "data_dirname" property of
    #BETSE rather than a subclass. One admirable means of circumventing this
    #might be as follows:
    #
    #* Redefine the data_yaml_dirname() property in terms of the
    #  betse_data_dirname() property (e.g., by substituting "self.data_dirname"
    #  for "self.betse_data_dirname").
    #* For both disambiguity and orthogonality, rename:
    #  * data_yaml_dirname() to betse_data_yaml_dirname().
    #  * sim_conf_default_filename() to betse_sim_conf_default_filename().
    #
    #Nice, eh? All existing semantics are preserved, including the capacity to
    #override these sane default implementations. Moreover, no additional giant
    #singleton need be preserved in memory merely to define a single property.

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
