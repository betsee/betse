#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Mandatory runtime dependency checking.

This module defines functions intended to be called by high-level interface
modules (e.g., `betse.cli.cli`) *before* attempting to import such dependencies.
'''

#FIXME: Implement setuptools-based version checking as well, for dependencies
#providing setuptools-specific egg metadata. stackoverflow offered the following
#useful example code:
#
#    from pkg_resources import (WorkingSet, DistributionNotFound)
#    working_set = WorkingSet()
#
#    # Printing all installed modules
#    print tuple(working_set)
#
#    # Detecting if module is installed
#    try:
#        dep = working_set.require('paramiko>=1.0')
#    except DistributionNotFound:
#        pass
#
#Given that, consider the following approach:
#
#    import metadata
#    import pkg_resources
#    from pkg_resources import (
#        DistributionNotFound, Environment, VersionConflict, WorkingSet)
#    working_set = WorkingSet()
#    environment = Environment(working_set.entries)
#    requirements = pkg_resources.parse_requirements(
#        metadata.REQUIREMENTS)
#    for requirement in requirements:
#        if not importlib.find_loader(requirement.project_name):
#            raise DistributionNotFound(
#                "Mandatory dependency {} not found.".format(requirement))
#            )
#
#        requirement_distribution = environment.best_match(
#            requirement, working_set)
#
#        # Oops, the "best" so far conflicts with a dependency.
#        if requirement_distribution and\
#           requirement_distribution not in requirements:
#               raise VersionConflict(
#                   "Mandatory dependency {} found but {} required.".format(
#                       requirement_distribution, requirement))
#
#The above *SHOULD* work, but assumes use of a new
#"metadata.install_requirements" module dict. Not terribly arduous, of course;
#merely shift such dict from its current use in "setup.py".

# ....................{ IMPORTS                            }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To raise human-readable exceptions on missing mandatory dependencies,
# the top-level of this module may import *ONLY* from packages guaranteed to
# exist at installation time (e.g., stock Python packages).
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

from betse import metadata
import importlib

# ....................{ EXCEPTIONS                         }....................
def die_unless_satisfied() -> None:
    '''
    Raise an exception if at least one mandatory `betse` dependency is unmet.

    This function *always* validates dependency existence. For dependencies with
    `setuptools`-specific metadata (e.g., `.egg-info/`-suffixed subdirectories
    of the `site-packages/` directory for the active Python 3 interpreter,
    typically created by `setuptools` at install time), this function *also*
    validates the version of such dependency.

    See Also
    ----------
    README.md
        Human-readable list of such dependencies.
    '''
    # If the setuptools-specific "pkg_resources" dependency is *NOT* importable,
    # fail *BEFORE* attempting to import such dependency below.
    die_unless_importable('pkg_resources', 'setuptools >= 7.0')

    # Import such dependency and all required classes in such dependency.
    from pkg_resources import (Environment, VersionConflict, WorkingSet)
    import pkg_resources

    # Set of all betse-specific dependencies in a manner usable below.
    requirements = pkg_resources.parse_requirements(metadata.REQUIREMENTS)

    # Set of all setuptools-installed Python packages available under the active
    # Python 3 interpreter.
    working_set = WorkingSet()

    # Helper for finding the former in the latter.
    environment = Environment(working_set.entries)

    # Validate each such dependency.
    for requirement in requirements:
        # If such dependency is *NOT* importable, fail.
        die_if_dependency_unloadable(requirement.project_name, str(requirement))

        # Else, such dependency is importable.
        #
        # Best match for such dependency under the active Python 3 interpreter.
        requirement_distribution = environment.best_match(
            requirement, working_set)

        # If the version of such match is *NOT* a version we require, fail.
        if requirement_distribution and\
           requirement_distribution not in requirements:
               raise VersionConflict(
                   'Mandatory dependency {} required but only {} found.'.format(
                       requirement, requirement_distribution))

def die_unless_importable(
    module_name: str,
    module_requirements: str) -> None:
    '''
    Raise an exception containing the passed human-readable module requirements
    if the module with the passed name is *not* importable.

    Such module is assumed to signify a mandatory `betse` dependency. Likewise,
    such requirements are assumed to be in `setuptools` format (e.g.,
    `numpy >= 1.9.0`).
    '''
    assert isinstance(module_name, str),\
        '"{}" not a string.'.format(module_name)
    assert isinstance(module_requirements, str),\
        '"{}" not a string.'.format(module_requirements)

    # If such dependency is *NOT* importable, raise an exception. Ideally, such
    # exception would be an instance of the setuptools-specific
    # "DistributionNotFound" class. Yet, as setuptools and hence such class is
    # *NOT* guaranteed to be importable, the conventional Python exception for
    # import errors is raised instead.
    if not importlib.find_loader(module_name):
        raise ImportError(
            'Mandatory dependency {} not found.'.format(module_requirements))

# --------------------( WASTELANDS                         )--------------------
        # if not importlib.find_loader(requirement.project_name):
        #     raise DistributionNotFound(
        #         'Mandatory dependency {} not found.'.format(requirement))
    # for module_name in {
    #     'scipy'
    # }:
    # For each such dependency, this function also attempts to validate such
    # dependency's version. Specifically, if such dependency both exists *and*,
    # an exception is raised if

    # Caveats
    # ----------
    # Dependency versions are *not* validated. This is subject to change
