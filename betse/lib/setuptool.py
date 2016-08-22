#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level support facilities for `pkg_resources`, a mandatory runtime
dependency simplifying inspection of application dependencies.
'''

# ....................{ IMPORTS                            }....................
import pkg_resources
from betse.exceptions import BetseDependencyException
from betse.util.type import sequences
from betse.util.type.types import type_check, MappingType, ModuleType
from collections import OrderedDict
from pkg_resources import (
    DistributionNotFound, Requirement, UnknownExtra, VersionConflict)

# ....................{ GLOBALS ~ dict                     }....................
SETUPTOOLS_TO_MODULE_NAME = {
    'Matplotlib': 'matplotlib',
    'Numpy': 'numpy',
    'Pillow': 'PIL',
    'PyYAML': 'yaml',
    'SciPy': 'numpy',
    'dill': 'dill',
    'networkx': 'networkx',
    'ptpython': 'ptpython',
    'pydot': 'pydot',
    'setuptools': 'setuptools',
    'six': 'six',
    'yamale': 'yamale',
}
'''
Dictionary mapping each relevant `setuptools`-specific project name (e.g.,
`PyYAML`) to the fully-qualified name of the corresponding top-level module or
package providing that project (e.g., `yaml`).

For consistency, the size of this dictionary should be greater than or equal to
the size of the `betse.metadata.DEPENDENCIES_RUNTIME_MANDATORY` unordered set.
'''

# ....................{ EXCEPTIONS                         }....................
@type_check
def die_unless_requirement_str(*requirement_strs: str) -> None:
    '''
    Raise an exception unless all passed `setuptools`-formatted requirements
    strings (e.g., `Numpy >= 1.8.0`) are satisfiable, implying the corresponding
    third-party packages to be both importable and of satisfactory version.

    See Also
    ----------
    :func:`is_requirement_str`
        Further details.
    '''

    # List of all high-level "Requirements" objects corresponding to the passed
    # low-level requirements strings.
    requirements = pkg_resources.parse_requirements(requirement_strs)

    # Validate these requirements.
    for requirement in requirements:
        die_unless_requirement(requirement)


@type_check
def die_unless_requirement(requirement: Requirement) -> None:
    '''
    Raise an exception unless the passed `setuptools`-specific requirement is
    satisfiable, implying the corresponding third-party package to be both
    importable and of satisfactory version.

    Parameters
    ----------
    requirement : Requirement
        `setuptools`-specific object encapsulating the required name and version
        of this package.
    '''

    # Human-readable exception to be raised below if any.
    exception_message = None

    # If setuptools raises a non-human-readable exception on validating this
    # requirement, convert that into a human-readable exception.
    try:
        # Validate this requirement via setuptools.
        pkg_resources.get_distribution(requirement)
    # If setuptools fails to find this requirement, this does *NOT* necessarily
    # imply this requirement to be unimportable as a package. See the
    # is_requirement() function for details.
    except DistributionNotFound:
        pass
    # If this requirement is found to be an insufficient version, a
    # non-human-readable exception resembling the following is raised:
    #
    #    pkg_resources.VersionConflict: (PyYAML 3.09 (/usr/lib64/python3.3/site-packages), Requirement.parse('PyYAML>=3.10'))
    except VersionConflict as version_conflict:
        exception_message = BetseDependencyException(
            'Dependency "{}" unsatisfied by installed dependency "{}".'.format(
                version_conflict.req, version_conflict.dist))
    #FIXME: Handle the "UnknownExtra" exception as well.

    # If a human-readable exception is to be raised, do so. While it would be
    # preferable to simply raise this exception in the exception handler above,
    # doing so encourages Python 3 to implicitly prepend this exception by the
    # non-human-readable exception raised above. (...Ugh!)
    if exception_message:
        raise exception_message

    # Avoid circular import dependencies.
    from betse.util.py import modules

    # Attempt to import this package.
    try:
        # print('Importing dependency: ' + module_name)
        package = import_requirement(requirement)
    # Convert this exception into a human-readable form.
    except ImportError:
        exception_message = BetseDependencyException(
            'Dependency "{}" not found.'.format(requirement.project_name))

    # If a human-readable exception is to be raised, do so.
    if exception_message:
        raise exception_message

    # Package version if any or "None" otherwise.
    package_version = modules.get_version_or_none(package)

    # If this version fails to satisfy this requirement, raise an exception.
    if package_version not in requirement:
        raise BetseDependencyException(
            'Dependency "{}" unsatisfied by installed version {}.'.format(
                requirement, package_version))

# ....................{ TESTERS                            }....................
@type_check
def is_requirement_str(*requirement_strs: str) -> bool:
    '''
    `True` only if all passed `setuptools`-formatted requirements strings (e.g.,
    `Numpy >= 1.8.0`) are satisfiable, implying the corresponding third-party
    packages to be both importable and of satisfactory version.

    Equivalently, this function returns `False` if at least one such
    requirement is unsatisfied. For importable unsatisfied dependencies with
    `setuptools`-specific metadata (e.g., `.egg-info/`-suffixed subdirectories
    of the `site-packages/` directory for the active Python 3 interpreter,
    typically created by `setuptools` at install time), this function
    additionally validates the versions of these dependencies to satisfy these
    requirements.
    '''

    # List of all high-level "Requirements" objects corresponding to these
    # low-level requirements strings.
    requirements = pkg_resources.parse_requirements(requirement_strs)

    # For each such requirement...
    for requirement in requirements:
        # If this requirement is unsatisfied, return "False".
        if not is_requirement(requirement):
            return False

    # Else, return "True".
    return True


@type_check
def is_requirement(requirement: Requirement) -> bool:
    '''
    `True` only if the passed `setuptools`-specific requirement is satisfiable,
    implying the corresponding third-party package to be both importable and of
    satisfactory version.

    Parameters
    ----------
    requirement : Requirement
        `setuptools`-specific object encapsulating this package's requisite name
        and version.

    Returns
    ----------
    bool
        `True` only if this requirement is satisfiable.
    '''

    # Attempt to...
    try:
        # Validate this requirement via setuptools.
        pkg_resources.get_distribution(requirement)

        # If no exception was raised, this requirement is satisfied.
        return True
    # If setuptools fails to find this requirement, this does *NOT* necessarily
    # imply this requirement to be unimportable as a package. Rather, this only
    # implies this requirement was *NOT* installed as a setuptools-managed egg.
    # This requirement is still installable and hence importable (e.g., by
    # manually copying this requirement's package into the "site-packages"
    # subdirectory of the top-level directory for this Python interpreter).
    # However, does this edge-case actually occur in reality? *YES.*
    # PyInstaller-frozen applications embed requirements without corresponding
    # setuptools-managed eggs. Hence, this edge-case *MUST* be handled.
    except DistributionNotFound:
        pass
    # If this requirement is found to be an insufficient version, fail.
    except (UnknownExtra, VersionConflict):
        return False

    # Avoid circular import dependencies.
    from betse.util.py import modules

    # If no setuptools-managed egg exists for this requirement, fallback to this
    # lower-level strategy:
    #
    # 1. Import this requirement's top-level package.
    # 2. Compare this package's "__version__" attribute (if any) with this
    #    requirement's required version (if any).
    #
    # Since this strategy is inherently less reliable than setuptools-based
    # dependency validation, the latter remains the default.
    # print('Validating dependency: ' + requirement.project_name)

    # Attempt to import this package.
    try:
        # print('Importing dependency: ' + package_name)
        package = import_requirement(requirement)
    # If this package is unimportable, fail.
    except ImportError:
        return False

    # Package version if any or "None" otherwise.
    package_version = modules.get_version_or_none(package)

    # Return "True" only if this version exists and satisfies this requirement.
    return package_version is not None and package_version in requirement

# ....................{ GETTERS                            }....................
@type_check
def get_requirement_str_metadata(*requirement_strs: str) -> OrderedDict:
    '''
    Ordered dictionary synopsizing the currently installed third-party packages
    corresponding to (but _not_ necessarily satisfying) the passed `setuptools`-
    formatted requirements strings (e.g., `Numpy >= 1.8.0`).

    For readability, these strings will be lexicographically presorted in
    ascending order.
    '''

    # Lexicographically sorted tuple of these strings.
    requirement_strs_sorted = sequences.sort_lexicographic_ascending(
        requirement_strs)

    # List of all high-level "Requirements" objects corresponding to these
    # low-level requirements strings.
    requirements = pkg_resources.parse_requirements(requirement_strs_sorted)

    # Ordered dictionary synopsizing these requirements
    metadata = OrderedDict()

    # For each such requirement...
    for requirement in requirements:
        # Name of this requirement.
        requirement_name = requirement.project_name

        # Version of this requirement.
        requirement_version = get_requirement_version_readable(requirement)

        # Append metadata describing this requirement to this dictionary.
        metadata[requirement_name + ' version'] = requirement_version

    # Return this dictionary.
    return metadata


@type_check
def get_requirement_version_readable(requirement: Requirement) -> str:
    '''
    Human-readable version string for the currently installed version of the
    third-party package corresponding to (but _not_ necessarily satisfying) the
    passed `setuptools`-specific requirement.

    This function is principally intended for use in printing package metadata
    in a non-critical manner and hence is guaranteed to _never_ raise fatal
    exceptions. If his package:

    * Is importable but fails to satisfy this requirement, a string describing
      this conflict is returned.
    * Is unimportable, the string `unimportable` is returned.
    * Has no `__version__` attribute, the string `unversioned` is returned.

    Parameters
    ----------
    requirement : Requirement
        `setuptools`-specific object encapsulating this package.

    Returns
    ----------
    str
        Version string for the currently installed version of this package if
        any or the string `unimportable` or `unversioned` (as detailed above).
    '''

    # Attempt to...
    try:
        # Query setuptools for the setuptools-specific "Distribution" object
        # describing the currently installed third-party package satisfying this
        # requirement if any.
        distribution = pkg_resources.get_distribution(requirement)

        # If no exception was raised, this requirement is satisfied. In this
        # case, return this package's version.
        return distribution.version
    # If setuptools fails to find this requirement, this does *NOT* necessarily
    # imply this requirement to be unimportable as a package. See the
    # is_requirement() function for details.
    except DistributionNotFound:
        pass
    # If this requirement is found to be an insufficient version, circumvent
    # this by returning this package's version as is.
    except VersionConflict as version_conflict:
        return '{} <fails to satisfy {}>'.format(
            version_conflict.dist.version,
            version_conflict.req)
    #FIXME: Handle the "UnknownExtra" exception as well.

    # Avoid circular import dependencies.
    from betse.util.py import modules

    # Attempt to import this package.
    try:
        # print('Importing dependency: ' + module_name)
        package = import_requirement(requirement)
    # If this package is unimportable, return an appropriate string.
    except ImportError:
        return 'unimportable'

    # Package version if any or "None" otherwise.
    package_version = modules.get_version_or_none(package)

    # If this version exists, return this version.
    if package_version is not None:
        return package_version
    # Else, return an appropriate string.
    else:
        return 'unversioned'

# ....................{ CONVERTERS                         }....................
@type_check
def convert_requirement_dict_to_strs(requirement_dict: MappingType) -> tuple:
    '''
    Convert the passed dictionary of `setuptools`-specific requirements strings
    into a tuple of such strings.

    This dictionary is assumed to map from the `setuptools`-specific project
    name of a third-party dependency (e.g., `NetworkX`) to the suffix of a
    `setuptools`-specific requirements string constraining this dependency
    (e.g., `>= 1.11`). Each element of the resulting tuple is a string of the
    form `{key} {value}`, converted from a key-value pair of this dictionary in
    arbitrary order.

    Parameters
    ----------
    requirement_dict : MappingType
        Dictionary of `setuptools`-specific requirements strings in the format
        described above.

    Returns
    ----------
    tuple
        Tuple of `setuptools`-specific requirements strings in the format
        described above.
    '''

    return convert_requirement_dict_keys_to_strs(
        requirement_dict, *requirement_dict.keys())


@type_check
def convert_requirement_dict_keys_to_strs(
    requirement_dict: MappingType, *requirement_names: str) -> tuple:
    '''
    Convert all key-value pairs of the passed dictionary of `setuptools`-
    specific requirements strings whose keys are the passed strings into a
    tuple of `setuptools`-specific requirements string.

    Parameters
    ----------
    requirement_dict : MappingType
        Dictionary of requirements strings.
    requirement_names : Tuple[str]
        Tuple of keys identifying the key-value pairs of this dictionary to
        convert.

    Returns
    ----------
    tuple
        Tuple of `setuptools`-specific requirements strings in the format
        described above.

    Raises
    ----------
    :exc:`BetseDependencyException`
        If the passed key is _not_ a key of this dictionary.

    See Also
    ----------
    :func:`convert_requirement_dict_to_strs`
        Further details on the format of this dictionary and resulting string.
    '''

    return tuple(
        convert_requirement_dict_key_to_str(requirement_dict, requirement_name)
        for requirement_name in requirement_names
    )


@type_check
def convert_requirement_dict_key_to_str(
    requirement_dict: MappingType, requirement_name: str) -> str:
    '''
    Convert the key-value pair of the passed dictionary of `setuptools`-
    specific requirements strings whose key is the passed string into a
    `setuptools`-specific requirements string.

    Parameters
    ----------
    requirement_dict : MappingType
        Dictionary of requirements strings.
    requirement_names : str
        Key identifying the key-value pairs of this dictionary to convert.

    Returns
    ----------
    str
        Requirements string converted from this key-value pair.

    Raises
    ----------
    :exc:`BetseDependencyException`
        If the passed key is _not_ a key of this dictionary.

    See Also
    ----------
    :func:`convert_requirement_dict_to_strs`
        Further details on the format of this dictionary and resulting string.
    '''

    # If this name is unrecognized, raise an exception.
    if requirement_name not in requirement_dict:
        raise BetseDependencyException(
            'Dependency "{}" unrecognized.'.format(requirement_name))

    # Convert this key-value pair into a requirements string.
    return '{} {}'.format(requirement_name, requirement_dict[requirement_name])

# ....................{ IMPORTERS                          }....................
@type_check
def import_requirement(requirement: Requirement) -> ModuleType:
    '''
    Import and return the top-level package object satisfying the passed
    `setuptools`-specific requirement.

    Parameters
    ----------
    requirement : Requirement
        `setuptools`-specific object encapsulating this package.

    Returns
    ----------
    ModuleType
        Top-level package object implementing this requirement.

    Raises
    ----------
    BetseDependencyException
        If this name is unrecognized (i.e., is _not_ a key of the
        `SETUPTOOLS_TO_MODULE_NAME` dictionary).
    ImportError
        If this package is unimportable.
    '''

    # Avoid circular import dependencies.
    from betse.util.py import modules

    # Name of this requirement.
    requirement_name = requirement.project_name

    # If this name is unrecognized, raise an exception.
    if requirement_name not in SETUPTOOLS_TO_MODULE_NAME:
        raise BetseDependencyException(
            'Requirement "{}" unrecognized.'.format(requirement_name))

    # Fully-qualified name of this requirement's package.
    package_name = SETUPTOOLS_TO_MODULE_NAME[requirement_name]

    # Import and return this package.
    return modules.import_module(package_name)
