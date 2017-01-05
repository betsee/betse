#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level support facilities for `pkg_resources`, a mandatory runtime
dependency simplifying inspection of application dependencies.
'''

# ....................{ IMPORTS                            }....................
import pkg_resources
from betse.exceptions import BetseLibException
from betse.util.type import iterables, modules
from betse.util.type.types import (
    type_check, GeneratorType, MappingType, ModuleType, NoneType, SequenceTypes)
from collections import OrderedDict
from pkg_resources import (
    Distribution,
    DistributionNotFound,
    Requirement,
    UnknownExtra,
    VersionConflict,
)

# ....................{ GLOBALS ~ dict                     }....................
SETUPTOOLS_TO_MODULE_NAME = {
    'Numpy': 'numpy',
    'Pillow': 'PIL',
    'Pympler': 'pympler',
    'PyYAML': 'yaml',
    'SciPy': 'numpy',
    'dill': 'dill',
    'matplotlib': 'matplotlib',
    'networkx': 'networkx',
    'numba': 'numba',
    'pprofile': 'pprofile',
    'ptpython': 'ptpython',
    'pydot': 'pydot',
    'pytest': 'pytest',
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

    Raises
    ----------
    BetseLibException
        If at least one such requirement is unsatisfiable.
    '''

    # List of all requirement objects parsed from these requirement strings.
    requirements = get_requirements(*requirement_strs)

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
        Object describing this package or module's required name and version.

    Raises
    ----------
    BetseLibException
        If this requirement is unsatisfiable.
    '''

    # Human-readable exception to be raised below if any.
    betse_exception = None

    try:
        # Object describing the currently installed version of the package or
        # module satisfying this requirement if any or "None" if this
        # requirement cannot be guaranteed to be unsatisfied.
        distribution = get_requirement_distribution_or_none(requirement)

        # If this requirement is satisfied, we're done here.
        if distribution is not None:
            return
        # Else, fallback to attempting to manually import this requirement.
    # If setuptools found only requirements of insufficient version, a
    # non-human-readable exception resembling the following is raised:
    #
    #    pkg_resources.VersionConflict: (PyYAML 3.09 (/usr/lib64/python3.3/site-packages), Requirement.parse('PyYAML>=3.10'))
    #
    # Detect this and raise a human-readable exception instead.
    except VersionConflict as version_conflict:
        betse_exception = BetseLibException(
            'Dependency "{}" unsatisfied by installed dependency "{}".'.format(
                version_conflict.req, version_conflict.dist))
    #FIXME: Handle the "UnknownExtra" exception as well.

    # If a human-readable exception is to be raised, do so. While it would be
    # preferable to simply raise this exception in the exception handler above,
    # doing so encourages Python 3 to implicitly prepend this exception by the
    # non-human-readable exception raised above. (...Ugh!)
    if betse_exception:
        raise betse_exception

    # Fully-qualified name of this requirement's package.
    package_name = get_requirement_module_name(requirement)

    # Attempt to manually import this package.
    try:
        # print('Importing dependency: ' + module_name)
        package = import_requirement(requirement)
    # If a standard import exception is raised...
    except ImportError as root_exception:
        # Low-level Python-specific exception message.
        root_exception_message = str(root_exception)

        # If this exception signifies the common case of a missing dependency,
        # avoid exposing this exception to end users. Doing so would convey no
        # meaningful metadata.
        if root_exception_message == (
            "No module named '{}'".format(package_name)):
            betse_exception = BetseLibException(
                'Dependency "{}" not found.'.format(requirement.project_name))
        # Else, this exception signifies an unexpected edge-case. For
        # debuggability, expose this exception to end users.
        else:
            raise BetseLibException(
                'Dependency "{}" unimportable.'.format(
                    requirement.project_name))
    # Else if any other exception is raised, expose this exception to end users.
    except Exception as root_exception:
        raise BetseLibException(
            'Dependency "{}" unimportable.'.format(requirement.project_name))

    # If a human-readable exception is to be raised, do so.
    if betse_exception:
        raise betse_exception

    # Package version if any or raise an exception otherwise.
    package_version = modules.get_version(package)

    # If this version fails to satisfy this requirement, raise an exception.
    if package_version not in requirement:
        raise BetseLibException(
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

    # List of all requirement objects parsed from these requirement strings.
    requirements = get_requirements(*requirement_strs)

    # If any such requirement is unsatisfied, fail.
    for requirement in requirements:
        if not is_requirement(requirement):
            return False

    # Else, succeed.
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
        Object describing this package or module's required name and version.

    Returns
    ----------
    bool
        `True` only if this requirement is satisfiable.
    '''

    try:
        # Object describing the currently installed version of the package or
        # module satisfying this requirement if any or "None" if this
        # requirement cannot be guaranteed to be unsatisfied.
        distribution = get_requirement_distribution_or_none(requirement)

        # If this requirement is satisfied, we're done here.
        if distribution is not None:
            return True
    # If setuptools found only requirements of insufficient version, fail.
    except (UnknownExtra, VersionConflict):
        return False

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
def get_requirements(*requirement_strs: str) -> GeneratorType:
    '''
    Generator of all high-level `setuptools`-specific
    :class:`pkg_resources.Requirement` objects parsed from the passed low-level
    requirement strings.

    Parameters
    ----------
    requirement_strs: Tuple[str]
        Tuple of all requirement strings to parse into requirement objects.

    Yields
    ----------
    Requirement
        For each passed requirement string, this generator yields a requirement
        object parsed from this string (_in the same order_)
    '''

    yield from pkg_resources.parse_requirements(requirement_strs)


@type_check
def get_requirement_distribution_or_none(
    requirement: Requirement) -> (Distribution, NoneType):
    '''
    Get a :class:`Distribution` instance describing the currently installed
    version of the top-level third-party package or module satisfying the passed
    `setuptools`-specific requirement if any _or_ raise an exception if this
    requirement is guaranteed to be unsatisfied (e.g., due to a version
    mismatch) _or_ `None` if this requirement cannot be guaranteed to be
    unsatisfied (e.g., due to this requirement being installed either without
    `setuptools` or with the `setuptools` subcommand `develop`).

    This high-level getter should _always_ be called in lieu of the low-level
    :func:`pkg_resources.get_distribution` function, which raises spurious
    exceptions in common non-erroneous edge cases (e.g., packages installed via
    the `setuptools` subcommand `develop`) and is thus unsafe for
    general-purpose use.

    Parameters
    ----------
    requirement : Requirement
        Object describing this package or module's required name and version.

    Returns
    ----------
    Distribution or NoneType
        Object describing the currently installed version of the package or
        module satisfying this requirement if any _or_ `None` otherwise.
        Specifically, `None` is returned in all of the following conditions --
        only one of which genuinely corresponds to an error:
        * This requirement is _not_ installed at all. (**Error.**)
        * This requirement was installed manually rather than with `setuptools`,
          in which case no such :class:`Distribution` exists. (**Non-error.**)
        * This requirement was installed with the `setuptools` subcommand
          `develop`, in which case a :class:`Distribution` technically exists
          but in a sufficiently inconsistent state that the low-level
          :func:`pkg_resources.get_distribution` function raises an exception on
          attempting to retrieve that object. (**Non-error.**)
        Since distinguishing the erroneous from non-erroneous cases exceeds the
        mandate of this getter, the caller is expected to do so.

    Raises
    ----------
    VersionConflict
        If the currently installed version of this package or module fails to
        satisfy this requirement's version constraints.
    '''

    # If setuptools finds this requirement, return its distribution as is.
    try:
        return pkg_resources.get_distribution(requirement)
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
        return None
    # If setuptools fails to find the distribution-packaged version of this
    # requirement (e.g., due to having been editably installed with "sudo
    # python3 setup.py develop"), this version may still be manually parseable
    # from this requirement's package. Since setuptools fails to raise an
    # exception whose type is unique to this error condition, the contents of
    # this exception are parsed to distinguish this expected error condition
    # from other unexpected error conditions. In the former case, a
    # non-human-readable exception resembling the following is raised:
    #
    #     ValueError: ("Missing 'Version:' header and/or PKG-INFO file", networkx [unknown version] (/home/leycec/py/networkx))
    except ValueError as version_missing:
        # If this exception was *NOT*...
        if not (
            # ...instantiated with two arguments
            len(version_missing.args) == 2 and
            # ...whose second argument is suffixed by a suffix indicating the
            # version of this distribution to have been ignored rather than
            # recorded during installation...
            str(version_missing.args[1]).endswith(' [unknown version]')
        ):
            # Then this exception indicates an unexpected and hence
            # non-ignorable error condition. Reraise this exception!
            raise

        # Else, this exception indicates an expected ignorable error condition.
        return None


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
    requirement_strs_sorted = iterables.sort_ascending(requirement_strs)

    # List of all requirement objects parsed from these requirement strings.
    requirements = get_requirements(*requirement_strs_sorted)

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
def get_requirement_module_name(requirement: Requirement) -> str:
    '''
    Fully-qualified name of the top-level module or package (e.g., :mod:`yaml`)
    providing the passed `setuptools`-specific requirement.

    Parameters
    ----------
    requirement : Requirement
        Object describing this package or module's required name and version.

    Returns
    ----------
    str
        Fully-qualified name of this module or package.

    Raises
    ----------
    BetseLibException
        If this name is unrecognized (i.e., is _not_ a key of the
        :data:`SETUPTOOLS_TO_MODULE_NAME` dictionary).
    '''

    # Name of this requirement.
    requirement_name = requirement.project_name

    # If this name is unrecognized, raise an exception.
    if requirement_name not in SETUPTOOLS_TO_MODULE_NAME:
        raise BetseLibException(
            'Requirement "{}" unrecognized.'.format(requirement_name))

    # Return the name of this requirement's module or package.
    return SETUPTOOLS_TO_MODULE_NAME[requirement_name]


@type_check
def get_requirement_version_readable(requirement: Requirement) -> str:
    '''
    Human-readable version string for the currently installed version of the
    third-party module or package corresponding to (but _not_ necessarily
    satisfying) the passed `setuptools`-specific requirement.

    This function is principally intended for use in printing package metadata
    in a non-critical manner and hence is guaranteed to _never_ raise fatal
    exceptions. If this package:

    * Is importable but fails to satisfy this requirement, a string describing
      this conflict is returned.
    * Is unimportable, the string `not installed` is returned.
    * Has no `__version__` attribute, the string `unknown` is returned.

    Parameters
    ----------
    requirement : Requirement
        Object describing this package or module's required name and version.

    Returns
    ----------
    str
        Version string for the currently installed version of this package if
        any or the string `not installed` or `unknown` (as detailed above).
    '''

    try:
        # Object describing the currently installed version of the package or
        # module satisfying this requirement if any or "None" if this
        # requirement cannot be guaranteed to be unsatisfied.
        distribution = get_requirement_distribution_or_none(requirement)

        # If this requirement is satisfied, return its version.
        if distribution is not None:
            return distribution.version
    # If setuptools found only requirements of insufficient version, return this
    # version regardless (with a suffix noting this to be the case).
    except VersionConflict as version_conflict:
        return '{} <fails to satisfy {}>'.format(
            version_conflict.dist.version,
            version_conflict.req)
    #FIXME: Handle the "UnknownExtra" exception as well.

    # Attempt to import this package.
    try:
        # print('Importing dependency: ' + module_name)
        package = import_requirement(requirement)
    # If this package is unimportable, return an appropriate string.
    except ImportError:
        return 'not installed'

    # Package version if any or "None" otherwise.
    package_version = modules.get_version_or_none(package)

    # If this version exists, return this version.
    if package_version is not None:
        return package_version
    # Else, return an appropriate string.
    else:
        return 'unknown'

# ....................{ CONVERTERS ~ tuple-to-dict         }....................
@type_check
def convert_requirements_tuple_to_dict(
    requirements_tuple: SequenceTypes) -> dict:
    '''
    Convert the passed tuple of `setuptools`-specific requirements strings into
    a dictionary of such strings.

    This tuple is assumed to contain strings of the form
    `{project_name} {project_specs}`, where:

    * `{project_name}` is the `setuptools`-specific project name of a
      third-party dependency (e.g., `NetworkX`).
    * `{project_specs}` is a `setuptools`-specific requirements string
      constraining this dependency (e.g., `>= 1.11`).

    Each key-value pair of the resulting dictionary maps from each such
    `{project_name}` to the corresponding `{project_specs}`.

    Parameters
    ----------
    requirements_tuple : SequenceTypes
        Tuple of `setuptools`-specific requirements strings in the format
        described above.

    Returns
    ----------
    dict
        Dictionary of `setuptools`-specific requirements strings in the format
        described above.
    '''

    # List of all requirement objects parsed from these requirement strings.
    requirements = get_requirements(*requirements_tuple)

    # Dictionary containing these requirements.
    requirements_dict = {}

    # For each such requirement...
    for requirement in requirements:
        # Comma-delimited string of all constraints of this requirement
        # iteratively reduced from the high-level "specs" sequence of 2-tuples
        # "(op, version)" of this requirement object into low-level strings.
        #
        # Technically, manually parsing each requirement string of the passed
        # tuple into the project name and specifications required below would
        # also be feasible. Since manual parsing is significantly more fragile
        # than deferring to the authoritative parsing already implemented by the
        # canonical "pkg_resources" module, however, the latter is preferred.
        requirement_specs_str = ','.join(
            '{} {}'.format(requirement_spec[0], requirement_spec[1])
            for requirement_spec in requirement.specs
        )

        # Add a new key-value pair to this dictionary mapping from the name of
        # this requirement to the comma-delimited string of all constraints of
        # this requirement.
        requirements_dict[requirement.project_name] = requirement_specs_str

    # Return this dictionary.
    return requirements_dict

# ....................{ CONVERTERS ~ dict-to-tuple         }....................
@type_check
def convert_requirements_dict_to_tuple(requirements_dict: MappingType) -> tuple:
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
    requirements_dict : MappingType
        Dictionary of `setuptools`-specific requirements strings in the format
        described above.

    Returns
    ----------
    tuple
        Tuple of `setuptools`-specific requirements strings in the format
        described above.
    '''

    return convert_requirements_dict_keys_to_tuple(
        requirements_dict, *requirements_dict.keys())


@type_check
def convert_requirements_dict_keys_to_tuple(
    requirements_dict: MappingType, *requirement_names: str) -> tuple:
    '''
    Convert all key-value pairs of the passed dictionary of `setuptools`-
    specific requirements strings whose keys are the passed strings into a
    tuple of `setuptools`-specific requirements strings.

    Parameters
    ----------
    requirements_dict : MappingType
        Dictionary of requirements strings.
    requirement_names : Tuple[str]
        Tuple of keys identifying the key-value pairs of this dictionary to
        convert.

    Returns
    ----------
    tuple
        Tuple of `setuptools`-specific requirements strings in the above format.

    Raises
    ----------
    :exc:`BetseLibException`
        If the passed key is _not_ a key of this dictionary.

    See Also
    ----------
    :func:`convert_requirements_dict_to_tuple`
        Further details on the format of this dictionary and resulting strings.
    '''

    return tuple(
        convert_requirements_dict_key_to_str(requirements_dict, requirement_name)
        for requirement_name in requirement_names
    )


@type_check
def convert_requirements_dict_key_to_str(
    requirements_dict: MappingType, requirement_name: str) -> str:
    '''
    Convert the key-value pair of the passed dictionary of `setuptools`-
    specific requirements strings whose key is the passed string into a
    `setuptools`-specific requirements string.

    Parameters
    ----------
    requirements_dict : MappingType
        Dictionary of requirements strings.
    requirement_names : str
        Key identifying the key-value pairs of this dictionary to convert.

    Returns
    ----------
    str
        Requirements string converted from this key-value pair.

    Raises
    ----------
    :exc:`BetseLibException`
        If the passed key is _not_ a key of this dictionary.

    See Also
    ----------
    :func:`convert_requirements_dict_to_tuple`
        Further details on the format of this dictionary and resulting string.
    '''

    # If this name is unrecognized, raise an exception.
    if requirement_name not in requirements_dict:
        raise BetseLibException(
            'Dependency "{}" unrecognized.'.format(requirement_name))

    # Convert this key-value pair into a requirements string.
    return '{} {}'.format(requirement_name, requirements_dict[requirement_name])

# ....................{ IMPORTERS                          }....................
@type_check
def import_requirement(requirement: Requirement) -> ModuleType:
    '''
    Import and return the top-level package object satisfying the passed
    `setuptools`-specific requirement.

    Parameters
    ----------
    requirement : Requirement
        Object describing this package or module's required name and version.

    Returns
    ----------
    ModuleType
        Top-level package object implementing this requirement.

    Raises
    ----------
    ImportError
        If this package is unimportable.
    '''

    # Fully-qualified name of this requirement's package.
    package_name = get_requirement_module_name(requirement)

    # Import and return this package.
    return modules.import_module(package_name)
