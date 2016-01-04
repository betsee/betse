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
    from betse.lib.matplotlib.mpl import mplconfig

    # Ensure that all mandatory dependencies exist *BEFORE* subsequent logic
    # (possibly) importing such dependencies.
    die_unless_satisfiable_all()

    # Configure such dependencies.
    mplconfig.init()

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
    from betse.util.python import modules

    # If the "pkg_resources" setuptools dependency is missing, raise an
    # exception *BEFORE* importing such dependency below.
    modules.die_unless(
        module_name = 'pkg_resources',
        exception_message = 'Mandatory dependency "pkg_resources" not found.'
    )

    # Validate such dependencies via "pkg_resources".
    from betse.lib import setuptool
    setuptool.die_unless_requirements_satisfiable_all()

# ....................{ GETTERS                            }....................
def get_metadata() -> OrderedDict:
    '''
    Get an ordered dictionary synopsizing all currently installed dependencies.
    '''
    # Imports deferred to their point of use, as documented above.
    from betse.lib import setuptool
    from betse.util.python import modules, pythons
    from betse.util.type import containers
    import pkg_resources

    # Dependency metadata to be collected and returned.
    dependency_metadata = OrderedDict()

    # If the active Python interpreter is frozen, query dependency versions
    # manually rather than via setuptools machinery. See the
    # setuptool.die_unless_requirement_satisfiable() function for related logic
    # and further commentary.
    if pythons.is_frozen():
    # if True:
        # List of the setuptools-specific project names of all BETSE
        # dependencies, lexicographically sorted for readability.
        project_names = containers.sort_as_lexicographic_ascending(
            setuptool.SETUPTOOLS_TO_MODULE_NAME.keys())

        # For each such name...
        for project_name in project_names:
            # Fully-qualified name of this project's root module or package.
            module_name =\
                setuptool.SETUPTOOLS_TO_MODULE_NAME[project_name]

	    # If this module is importable and hence frozen with this
	    # executable, this module is a describable dependency.
            if modules.is_module(module_name):
                # Version specifier provided by that module or package.
                module_version = modules.get_version(module_name)

                # Append metadata describing such dependency.
                dependency_metadata[project_name + ' version'] = module_version
    # Else, the active Python interpreter is *NOT* frozen. In such case, query
    # dependency versions via the more reliable setuptools machinery.
    else:
        # List of all BETSE dependencies as setuptools-specific requirements,
        # lexicographically sorted for readability.
        requirements = pkg_resources.parse_requirements(
            containers.sort_as_lexicographic_ascending(
                metadata.DEPENDENCIES_RUNTIME))

        # For each such dependency...
        for requirement in requirements:
            # Setuptools distribution describing such dependency. Since the
            # previously called dependencies.init() function presumably
            # succeeded, such distribution is guaranteed to exist.
            distribution = pkg_resources.get_distribution(requirement)

            # Append metadata describing such dependency.
            dependency_metadata[distribution.project_name + ' version'] =\
                distribution.version

    # Get such dictionary.
    return dependency_metadata

# --------------------( WASTELANDS                         )--------------------
        # Human-readable exception to be raised below if any.
        # exception = None

        # List of the fully-qualified names of the top-level modules and
        # packages specific to all BETSE dependencies, lexicographically sorted
        # for readability.
        # module_names = containers.sort_as_lexicographic_ascending(
        #     metadata.DEPENDENCY_TO_MODULE_NAME.values()

        # If the current dependency is setuptools *AND* the active Python
        # interpreter is frozen, avoid attempting to validate this dependency.
        # Since the BETSE codebase imports the setuptools-provided
        # "pkg_resources" package but *NOT* the separate "setuptools" package,
        # BETSE is usually frozen with the former rather than the latter. Since
        # the existence of the former has already been validated above, assume
        # for simplicity that the currently frozen version of "pkg_resources"
        # satisfies this dependency.
        # print('Validating dependency: ' + requirement.project_name)
        # if requirement.project_name == 'setuptools' and pythons.is_frozen():
        #     # print('Skipped!')
        #     continue

#FUXME: Ugh. This doesn't appear to work as expected, as trivially verified by
#the fact that the prior "yaml >= 3.10" specification appeared to work when that
#should have read "pyyaml >= 3.10". Even if this function is (somehow)
#mystically correct, however, it's still overblown. Given what we now know of
#setuptools, this is reducible to the following algorithm:
#
#* Iterate through the list of "betse.metadata.DEPENDENCIES_RUNTIME" setuptools
#  specifications.
#* For each such specification, pass that to the
#  pkg_resources.get_distribution() function.
#* That returns the setuptools distribution for that dependency. Given that,
#  simply perform the following test and raise an exception as required:
#
#    if dependency_distro not in pkg_resources.parse_version(depndency_spec):
#        raise VersionConflict('yadda yadda')
#
#Done. Awesome it is.

#FUXME: Terrible. This is implementable with pure-setuptools-based logic, which
#would have the distinct advantage of being implicitly synchronized with the
#list of required dependencies in "betse.metadata" *WITHOUT* explicitly
#requiring that this function ever be manually changed. To do so:
#
#* Iterate through the list of "betse.metadata.DEPENDENCIES_RUNTIME" setuptools
#  specifications.
#* For each such specification, pass that to the
#  pkg_resources.get_distribution() function.
#* That returns the setuptools distribution for that dependency. Given that,
#  simply append the following tuple to the returned OrderedDict:
#    (distro.project_name + ' version', distro.version)
#
#Done. Awesome it is.

        # try:
        # except VersionConflict as version_conflict:
    # # Imports deferred to their point of use, as documented above.
    # import matplotlib, numpy, scipy
    #
    # # Get such dictionary.
    # return OrderedDict((
    #     ('matplotlib version', matplotlib.__version__),
    #     ('numpy version', numpy.__version__),
    #     ('scipy version', scipy.__version__),
    # ))

# def die_unless_satisfiable() -> None:
#     '''
#     Raise an exception unless mandatory runtime dependencies of `betse` are
#     **satisfiable** (i.e., importable and of a satisfactory version).
#
#     Equivalently, an exception is raised if at least one such dependency is
#     unsatisfied.
#
#     Specifically, this function unconditionally validates the existence of all
#     such dependencies. For such dependencies with `setuptools`-specific metadata
#     (e.g., `.egg-info/`-suffixed subdirectories of the `site-packages/`
#     directory for the active Python 3 interpreter, typically created by
#     `setuptools` at install time), this function additionally validates the
#     versions of such dependencies to satisfy `betse` requirements.
#     '''
#     # Template for exception messages raised on missing dependencies.
#     exception_template = 'Mandatory dependency "{}" not found.'
#
#     # If the setuptools-specific "pkg_resources" dependency is missing, fail
#     # *BEFORE* attempting to import such dependency below.
#     modules.die_unless(
#         module_name = 'pkg_resources',
#         exception_message = exception_template.format('pkg_resources')
#     )
#
#     # Import such dependency and all required classes in such dependency.
#     from pkg_resources import (Environment, VersionConflict, WorkingSet)
#     import pkg_resources
#
#     # Set of all BETSE-specific dependencies as instances of the "Requirements"
#     # setuptools class.
#     requirements = pkg_resources.parse_requirements(
#         metadata.DEPENDENCIES_RUNTIME)
#
#     # Set of all setuptools-installed Python packages available under the active
#     # Python 3 interpreter.
#     # working_set = WorkingSet()
#
#     # Helper for finding the former in the latter.
#     # environment = Environment(pkg_resources.working_set.entries)
#
#     # Validate each such dependency.
#     for requirement_required in requirements:
#         # Name of the top-level importable module provided by such project.
#         module_name = requirement_required.project_name
#
#         # If such dependency is missing, fail.
#         if not modules.is_module(module_name):
#             # If such dependency is "setuptools", reduce such dependency to
#             # merely "pkg_resources" and try again. "pkg_resources" is a single
#             # module installed with (but inexplicably outside of the package
#             # tree of) setuptools. BETSE requires setuptools and hence
#             # "pkg_resources" at install time but *ONLY* the latter at runtime.
#             if module_name == 'setuptools':
#                 module_name = 'pkg_resources'
#
#             # Try again.
#             modules.die_unless(
#                 module_name, exception_template.format(
#                     str(requirement_required)))
#
#         # Else, such dependency exists.
#         #
#         # Best match for such dependency under the active Python 3 interpreter
#         # if any or None otherwise.
#         requirement_provided = environment.best_match(
#             requirement_required, working_set)
#
#         # If a match was found with version conflicting with that required,
#         # fail.
#         if requirement_provided is not None and\
#            requirement_provided not in requirement_required:
#                raise VersionConflict(
#                    'Mandatory dependency {} required but only {} found.'.format(
#                        requirement_required, requirement_provided))

    # Import such dependency and all required classes in such dependency.
    # from pkg_resources import ,
#To do this function unconditionally validates the existence of all
    # such dependencies.
#FUXME: Refactor as follows:
#
#* Create a new subpackage "betse.util.dependency".
#* Split this module into two modules:
#  * "betse.util.dependency.dependencies", providing the exception handling
#    function.
#  * "betse.util.dependency.matplotlib", providing the matplotlib-specific
#    functions.
#
#Hence, this module survives, albeit in *EXTREMELY* limited form. That's fine,
#for now. (Better minimal than overkill, we should think.)

    # If the current operating system is Apple OS X, prefer the "CocoaAgg"
    # backend to the "MacOSX" backend. The former leverages the cross-platform
    # C++ library AGG (Anti-grain Geometry) and hence tends to be better
    # supported; the latter does not.
        #FUXME: Extract into a new packages.is_package_PyObjC() function.
        # If PyObjC is installed, enable the "CocoaAgg" backend, which
        # internally requires such dependency.
        #if modules.is_module('PyObjCTools'):
#           loggers.log_warning(
#               'Optional dependency "PyObjC" not found. '
#               'Falling back from matplotlib backend "CocoaAgg" to "MacOSX".'
#           )
#           matplotlib.use('macosx')

        #FUXME: Alternately, perhaps we want to redirect
            # exception_template.format(metadata.DEPENDENCY_SETUPTOOLS)
    # List of setuptools requirements strings signifying all safely testable
    # mandatory dependencies. Sadly, the following mandatory dependencies are
    # *NOT* safely testable:
    #
    # * PySide. Under numerous Linux distributions, PySide is installed with
    #   cmake rather than setuptools
    # 'pyside >= 1.2.0',

#FUXME: Implement setuptools-based version checking as well, for dependencies
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
