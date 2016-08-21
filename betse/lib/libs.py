#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level **dependency** (i.e., both mandatory and optional third-party Python
packages imported at runtime) facilities.

This module defines functions intended to be called by high-level interface
modules (e.g., `betse.cli.clicli`) _before_ attempting to import dependencies.
'''

# ....................{ IMPORTS                            }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To raise human-readable exceptions on missing mandatory dependencies,
# the top-level of this module may import *ONLY* from packages guaranteed to
# exist at installation time (i.e., stock Python and BETSE packages).
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

from betse import metadata
from collections import OrderedDict

# ....................{ EXCEPTIONS                         }....................
def die_unless_satisfied_runtime_mandatory_all() -> None:
    '''
    Raise an exception unless all mandatory runtime dependencies of BETSE are
    **satisfiable** (i.e., importable and of a satisfactory version).

    Equivalently, this function raises an exception if at least one such
    dependency is unsatisfied. For importable unsatisfied dependencies with
    `setuptools`-specific metadata (e.g., `.egg-info/`-suffixed subdirectories
    of the `site-packages/` directory for the active Python 3 interpreter,
    typically created by `setuptools` at install time), this function
    additionally validates the versions of such dependencies to satisfy `betse`
    requirements.
    '''

    # Avoid circular import dependencies.
    from betse.util.py import modules

    # If the "pkg_resources" setuptools dependency is missing, raise an
    # exception *BEFORE* importing this dependency below.
    modules.die_unless_module(
        module_name='pkg_resources',
        exception_message='Mandatory dependency "pkg_resources" not found.',
    )

    # Validate these dependencies via "pkg_resources". Defer the importation of
    # this submodule until *AFTER* validating "pkg_resources" to exist above.
    from betse.lib import setuptool
    setuptool.die_unless_requirement_str(
        *metadata.DEPENDENCIES_RUNTIME_MANDATORY)

# ....................{ INITIALIZERS                       }....................
def init() -> None:
    '''
    Initialize all mandatory runtime dependencies of BETSE.

    Specifically, this function:

    * Reconfigures matplotlib with sane defaults specific to the current
      platform.
    '''

    # Avoid circular import dependencies.
    from betse.lib.matplotlib.matplotlibs import mpl_config
    from betse.lib.yaml import yamls

    # Configure these dependencies.
    mpl_config.init()
    yamls.init()

# ....................{ GETTERS                            }....................
def get_runtime_mandatory_metadata() -> OrderedDict:
    '''
    Ordered dictionary describing all currently installed third-party
    dependencies required by core functionality at runtime.
    '''

    # Avoid circular import dependencies.
    from betse.lib import setuptool

    # Return this metadata.
    return setuptool.get_requirement_str_metadata(
        *metadata.DEPENDENCIES_RUNTIME_MANDATORY)


def get_runtime_optional_metadata() -> OrderedDict:
    '''
    Ordered dictionary describing all currently installed third-party
    dependencies required by optional functionality at runtime.
    '''

    # Avoid circular import dependencies.
    from betse.lib import setuptool

    # Return this metadata.
    return setuptool.get_requirement_str_metadata(
        *metadata.DEPENDENCIES_RUNTIME_OPTIONAL)


def get_testing_mandatory_metadata() -> OrderedDict:
    '''
    Ordered dictionary describing all currently installed third-party
    dependencies required by this application's test suite.
    '''

    # Avoid circular import dependencies.
    from betse.lib import setuptool

    # Return this metadata.
    return setuptool.get_requirement_str_metadata(
        *metadata.DEPENDENCIES_TESTING_MANDATORY)
