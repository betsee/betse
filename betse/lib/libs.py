#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level **dependency** (i.e., both mandatory and optional Python packages
imported at runtime) facilities.

This module provides functions intended to be called by high-level interface
modules (e.g., `betse.cli.cli`) *before* attempting to import such dependencies.
'''

# ....................{ IMPORTS                            }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To raise human-readable exceptions on missing mandatory dependencies,
# the top-level of this module may import *ONLY* from packages guaranteed to
# exist at installation time (i.e., stock Python and BETSE packages).
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

from collections import OrderedDict
from betse import metadata

# ....................{ INITIALIZERS                       }....................
def init() -> None:
    '''
    Initialize all mandatory runtime dependencies of `betse`.

    This function (in order):

    . Raises an exception unless all dependencies are currently satisfiable.
    . Reconfigure `matplotlib` with sane defaults specific to the current
      platform.
    '''

    # Avoid circular import dependencies.
    from betse.lib.matplotlib.matplotlibs import mpl_config

    # Ensure that all mandatory dependencies exist *BEFORE* subsequent logic
    # (possibly) importing such dependencies.
    die_unless_satisfiable_all()

    # Configure such dependencies.
    mpl_config.init()

# ....................{ EXCEPTIONS                         }....................
def die_unless_satisfiable_all() -> None:
    '''
    Raise an exception unless all mandatory runtime dependencies of `betse` are
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
    # exception *BEFORE* importing such dependency below.
    modules.die_unless_module(
        module_name='pkg_resources',
        exception_message='Mandatory dependency "pkg_resources" not found.',
    )

    # Validate these dependencies via "pkg_resources".
    from betse.lib import setuptool
    setuptool.die_unless_requirements_satisfiable_all()

# ....................{ GETTERS                            }....................
def get_metadata() -> OrderedDict:
    '''
    Get an ordered dictionary synopsizing all currently installed dependencies.
    '''

    # Imports deferred to their point of use, as documented above.
    import pkg_resources
    from betse.lib import setuptool
    from betse.util.py import modules, pys
    from betse.util.type import sequences

    # Dependency metadata to be collected and returned.
    dependency_metadata = OrderedDict()

    # If the active Python interpreter is frozen, query dependency versions
    # manually rather than via setuptools machinery. See the
    # setuptool.die_unless_requirement_satisfiable() function for related logic
    # and further commentary.
    if pys.is_frozen():
    # if True:
        # List of the setuptools-specific project names of all BETSE
        # dependencies, lexicographically sorted for readability.
        project_names = sequences.sort_lexicographic_ascending(
            setuptool.SETUPTOOLS_TO_MODULE_NAME.keys())

        # For each such name...
        for project_name in project_names:
            # Fully-qualified name of this project's root module or package.
            module_name = (
                setuptool.SETUPTOOLS_TO_MODULE_NAME[project_name])

            # If this module is importable and hence frozen with this
            # executable, this module is a describable dependency.
            if modules.is_module(module_name):
                # Version specifier provided by that module or package.
                module_version = modules.get_version(module_name)

                # Append metadata describing this dependency.
                dependency_metadata[project_name + ' version'] = module_version
    # Else, the active Python interpreter is *NOT* frozen. In such case, query
    # dependency versions via the more reliable setuptools machinery.
    else:
        # List of all BETSE dependencies as setuptools-specific requirements,
        # lexicographically sorted for readability.
        requirements = pkg_resources.parse_requirements(
            sequences.sort_lexicographic_ascending(
                metadata.DEPENDENCIES_RUNTIME))

        # For each such dependency...
        for requirement in requirements:
            # Setuptools distribution describing such dependency. Since the
            # previously called dependencies.init() function presumably
            # succeeded, this distribution is guaranteed to exist.
            distribution = pkg_resources.get_distribution(requirement)

            # Append metadata describing this dependency.
            dependency_metadata[distribution.project_name + ' version'] = (
                distribution.version)

    # Return this dictionary.
    return dependency_metadata
