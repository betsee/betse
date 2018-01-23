#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level support facilities for :mod:`pkg_resources`, a mandatory runtime
dependency simplifying inspection of application dependencies.
'''

# ....................{ IMPORTS                            }....................
import pkg_resources
from betse.exceptions import BetseLibException
from betse.util.io.log import logs
from betse.util.type import iterables, modules
from betse.util.type.modules import DISTUTILS_PROJECT_NAME_TO_MODULE_NAME
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

# ....................{ EXCEPTIONS                         }....................
@type_check
def die_unless_requirements_dict(requirements_dict: MappingType) -> None:
    '''
    Raise an exception unless each :mod:`setuptools`-formatted requirements
    string produced by concatenating each key and corresponding value of the
    passed dictionary (e.g., converting key ``Numpy`` and value ``>= 1.8.0``
    into requirements string ``Numpy >= 1.8.0``) are satisfiable, implying these
    third-party packages to be both importable and of satisfactory version.

    Parameters
    ----------
    requirements_dict : MappingType
        Dictionary of requirements strings to validate.

    Raises
    ----------
    BetseLibException
        If at least one such requirement is unsatisfiable.
    '''

    # Tuple of requirements strings converted by concatenating each key and
    # value of this dictionary.
    requirement_strs = convert_requirements_dict_to_tuple(requirements_dict)

    # Validate these requirements strings.
    die_unless_requirement_str(*requirement_strs)


@type_check
def die_unless_requirements_dict_keys(
    requirements_dict: MappingType, *requirement_names: str) -> None:
    '''
    Raise an exception unless each :mod:`setuptools`-formatted requirements
    string produced by concatenating each passed key and corresponding value of
    the passed dictionary (e.g., converting key ``Numpy`` and value ``>= 1.8.0``
    into requirements string ``Numpy >= 1.8.0``) are satisfiable, implying these
    third-party packages to be both importable and of satisfactory version.

    Parameters
    ----------
    requirements_dict : MappingType
        Dictionary of requirements strings to validate a subset of.
    requirement_names : Tuple[str]
        Tuple of keys identifying the key-value pairs of this dictionary to
        validate.

    Raises
    ----------
    BetseLibException
        If at least one such requirement is unsatisfiable.
    '''

    # Tuple of requirements strings converted by concatenating each such key and
    # value of this dictionary.
    requirement_strs = convert_requirements_dict_keys_to_tuple(
        requirements_dict, *requirement_names)

    # Validate these requirements strings.
    die_unless_requirement_str(*requirement_strs)


@type_check
def die_unless_requirement_str(*requirement_strs: str) -> None:
    '''
    Raise an exception unless all passed :mod:`setuptools`-formatted
    requirements strings (e.g., ``Numpy >= 1.8.0``) are satisfiable, implying
    the corresponding third-party packages to be both importable and of
    satisfactory version.

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
    Raise an exception unless the passed :mod:`setuptools`-specific requirement
    is satisfiable, implying the corresponding third-party package to be both
    importable and of satisfactory version.

    Parameters
    ----------
    requirement : Requirement
        Object describing this module or package's required name and version.

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

    # If raising a human-readable exception, do so. While raising this exception
    # in the exception handler above would be preferable, doing so encourages
    # Python 3 to implicitly prepend this exception by the non-human-readable
    # exception raised above. This convoluted logic circumvents that. (...Ugh!)
    if betse_exception:
        raise betse_exception

    # Fully-qualified name of this requirement's package.
    package_name = DISTUTILS_PROJECT_NAME_TO_MODULE_NAME[
        requirement.project_name]

    # Attempt to manually import this requirement's package.
    try:
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

    # If this requirement is unversioned, return immediately for safety. While
    # the __contains__() methods of unversioned requirements implicitly called
    # by the "in" operator below do correctly return True for all possible
    # versions, this package may *NOT* declare the optional "__version__"
    # attribute; in that case, calling the modules.get_version() function below
    # would raise a fatal exception. Avoid this by short-circuiting immediately.
    if not _is_requirement_versioned(requirement):
        return
    # Else, this requirement is versioned.

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
    ``True`` only if all passed :mod:`setuptools`-formatted requirement strings
    (e.g., `Numpy >= 1.8.0`) are satisfiable, implying the corresponding
    third-party packages to be importable *and* of satisfactory versions.

    Parameters
    ----------
    requirement_strs : tuple[str]
        Tuple of all requirement strings to test.

    Returns
    ----------
    bool
        ``True`` only if these requirements are all satisfiable.
    '''

    # List of all requirement objects parsed from these requirement strings.
    requirements = get_requirements(*requirement_strs)

    # Return "True" only if these requirements are all satisfiable.
    return all(
        is_requirement(requirement)
        for requirement in requirements
    )


@type_check
def is_requirement(requirement: Requirement) -> bool:
    '''
    ``True`` only if the passed :mod:`setuptools`-specific **requirement**
    (i.e., object describing this module or package's required name and version)
    is satisfiable, implying the corresponding third-party package to be both
    importable and of satisfactory version.

    Parameters
    ----------
    requirement : Requirement
        Requirement to test.

    Returns
    ----------
    bool
        ``True`` only if this requirement is satisfiable.
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

    # Attempt to manually import this requirement's package.
    try:
        package = import_requirement(requirement)
    # If this package is unimportable, fail.
    except ImportError:
        return False

    # If this requirement is unversioned, all possible versions of this package
    # satisfy this requirement, in which case this requirement is satisfied.
    if not _is_requirement_versioned(requirement):
        return True
    # Else, this requirement is versioned.

    # Package version if any or "None" otherwise.
    package_version = modules.get_version_or_none(package)

    # Return "True" only if this version exists and satisfies this requirement.
    return package_version is not None and package_version in requirement

# ....................{ TESTERS ~ private                  }....................
@type_check
def _is_requirement_versioned(requirement: Requirement) -> bool:
    '''
    ``True`` only if the passed :mod:`setuptools`-specific requirement is
    **versioned** (i.e., constrained to require only a subset of all available
    versions of this requirement's distribution).

    This

    Parameters
    ----------
    requirement : Requirement
        Object describing this module or package's required name and version.

    Returns
    ----------
    bool
        ``True`` only if this requirement is versioned.
    '''

    # The "Requirement.specifier" instance variable of type "SpecifierSet"
    # encapsulates the set of all versions required by this requirement. Note:
    #
    # * The SpecifierSet.__init__() method accepts a string defining this set
    #   (e.g., ">= 3.1415"), which is the empty string when no specific versions
    #   are required, in which case this set is empty.
    # * The SpecifierSet.__eq__() method permits such sets to compared against
    #   such strings.
    #
    # Together, these two methods imply that a requirement is versioned if and
    # only if this requirement's specifier is *NOT* the empty string.
    return requirement.specifier != ''

# ....................{ GETTERS                            }....................
@type_check
def get_requirements(*requirement_strs: str) -> GeneratorType:
    '''
    Generator of all high-level :mod:`setuptools`-specific
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
    :class:`Distribution` instance describing the currently installed version of
    the top-level third-party module or package satisfying the passed
    :mod:`setuptools`-specific requirement if any *or* raise an exception if
    this requirement is guaranteed to be unsatisfied (e.g., due to a version
    mismatch) *or* ``None`` if this requirement cannot be guaranteed to be
    unsatisfied (e.g., due to this requirement being installed either without
    :mod:`setuptools` or with the :mod:`setuptools` subcommand ``develop``).

    This high-level getter should *always* be called in lieu of the low-level
    :func:`pkg_resources.get_distribution` function, which raises spurious
    exceptions in common non-erroneous edge cases (e.g., packages installed via
    the :mod:`setuptools` subcommand ``develop``) and is thus unsafe for
    general-purpose use.

    Parameters
    ----------
    requirement : Requirement
        Object describing this module or package's required name and version.

    Returns
    ----------
    Distribution or NoneType
        Object describing the currently installed version of the package or
        module satisfying this requirement if any *or* ``None`` otherwise.
        Specifically, ``None`` is returned in all of the following conditions --
        only one of which genuinely corresponds to an error:
        * This requirement is *not* installed at all. (**Error.**)
        * This requirement was installed manually rather than with
          :mod:`setuptools`, in which case no such :class:`Distribution` exists.
          (**Non-error.**)
        * This requirement was installed with the :mod:`setuptools` subcommand
          ``develop``, in which case a :class:`Distribution` technically exists
          but in a sufficiently inconsistent state that the low-level
          :func:`pkg_resources.get_distribution` function raises an exception on
          attempting to retrieve that object. (**Non-error.**)
        Since distinguishing the erroneous from non-erroneous cases exceeds the
        mandate of this getter, the caller is expected to do so.

    Raises
    ----------
    VersionConflict
        If the currently installed version of this module or package fails to
        satisfy this requirement's version constraints.
    '''

    # If setuptools finds this requirement, return its distribution as is.
    try:
        # logs.log_debug(
        #     'Retrieving requirement "%r" distribution...', requirement)
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
        # logs.log_debug(
        #     'Requirement "%r" distribution not found.', requirement)
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
        # If this exception was...
        if (
            # ...instantiated with two arguments...
            len(version_missing.args) == 2 and
            # ...whose second argument is suffixed by a string indicating the
            # version of this distribution to have been ignored rather than
            # recorded during installation...
            str(version_missing.args[1]).endswith(' [unknown version]')
        # ...this exception indicates an expected ignorable error condition.
        # Silently ignore this edge case.
        ):
            # logs.log_debug(
            #     'Requirement "%r" distribution version not found.', requirement)
            return None

        # Else, this exception indicates an unexpected and hence
        # non-ignorable error condition. Reraise this exception!
        raise

# ....................{ GETTERS ~ requirement : metadata   }....................
@type_check
def get_requirement_metadata(requirement: Requirement) -> str:
    '''
    Human-readable string describing the currently installed third-party
    module or package corresponding to (but *not* necessarily satisfying) the
    passed :mod:`setuptools`-specific requirement.

    This function is principally intended for use in printing package metadata
    in a non-critical manner and hence is guaranteed to *never* raise fatal
    exceptions. If this module or package:

    * Is unimportable, the string ``not installed`` is returned.
    * Is importable and...
      * Fails to satisfy this requirement, a string describing this conflict is
        returned.
      * Satisfies this requirement, the string ``{version} <{pathname}>`` is
        returned, where:
        * ``{pathname}`` is the absolute path of this module or package.
        * ``{version}`` is either:
          * If this module or package exports a version (e.g., via the
            ``__version__`` attribute), this version as a `.`-delimited string.
          * Else, the string ``unknown version``.

    Parameters
    ----------
    requirement : Requirement
        Object describing this module or package's required name and version.

    Returns
    ----------
    str
        Human-readable string describing this requirement.
    '''

    # Distribution satisfying this requirement if any or "None" otherwise.
    distribution = None

    try:
        # Object describing the currently installed version of the package or
        # module satisfying this requirement if any or "None" if this
        # requirement cannot be guaranteed to be unsatisfied.
        distribution = get_requirement_distribution_or_none(requirement)
    # If setuptools found only requirements of insufficient version, return this
    # version regardless (with a suffix noting this to be the case).
    except VersionConflict as version_conflict:
        return '{} [fails to satisfy {}]'.format(
            version_conflict.dist.version, version_conflict.req)
    #FIXME: Handle the "UnknownExtra" exception as well.

    # Attempt to manually import this module or package.
    try:
        package = import_requirement(requirement)
    # If this package is unimportable, return a human-readable string.
    except ImportError:
        return 'not installed'

    # Pathname and version of this module or package.
    package_pathname = modules.get_filename(package)
    package_version = None

    # If this requirement is satisfied, reuse the version provided by this
    # requirement's high-level distribution -- which is guaranteed to be at
    # least as precise as version metadata provided by this requirement's
    # low-level module or package (if any).
    #
    # Unfortunately, although the "Distribution" class does provide a public
    # "location" instance variable, this variable's value is typically the
    # absolute path of the parent directory containing this module or package
    # rather than the absolute path of the latter. The former is overly
    # ambiguous and hence useless for our purposes. We have no recourse but to
    # manually import this module or package and inspect its attributes. (*UGH*)
    if distribution is not None:
        package_version = distribution.version
    # Else, this requirement is unsatisfied. Fallback to the version metadata
    # provided by this requirement's low-level module or package.
    else:
        package_version = modules.get_version_or_none(package)

        # If no such version is provided, default to a human-readable string.
        if package_version is None:
            package_version = 'unknown version'

    # Return the expected string in the event of success.
    return '{} <{}>'.format(package_version, package_pathname)

# ....................{ GETTERS ~ requirements : metadata  }....................
@type_check
def get_requirements_dict_metadata(
    requirements_dict: MappingType) -> OrderedDict:
    '''
    Ordered dictionary synopsizing the currently installed third-party packages
    corresponding to (but *not* necessarily satisfying) the
    :mod:`setuptools`-formatted requirements strings produced by concatenating
    each key and value of the passed dictionary (e.g., converting key ``Numpy``
    and value ``>= 1.8.0`` into requirements string ``Numpy >= 1.8.0``).

    Parameters
    ----------
    requirements_dict : MappingType
        Dictionary of requirements strings to retrieve metadata for.

    Returns
    ----------
    OrderedDict
        Ordered dictionary lexicographically presorted on keys in ascending
        order such that each:
        * Key is a string of the form ``{{requirement_name}} version``.
        * Value is either:
          * The currently installed version specifier of this requirement
            if this requirement is both installed and exports an
            inspectable version at runtime (e.g., an ``__version__``
            package attribute).
          * The string ``not installed`` if this requirement is not found.
          * The string ``unknown version`` if this requirement exports no
            inspectable version at runtime.
    '''

    # Tuple of requirements strings converted by concatenating each key and
    # value of this dictionary.
    requirement_strs = convert_requirements_dict_to_tuple(requirements_dict)

    # Return this metadata.
    return get_requirement_str_metadata(*requirement_strs)


@type_check
def get_requirement_str_metadata(*requirement_strs: str) -> OrderedDict:
    '''
    Ordered dictionary synopsizing the currently installed third-party packages
    corresponding to (but *not* necessarily satisfying) the passed
    :mod:`setuptools`-formatted requirements strings (e.g., ``Numpy >= 1.8.0``).

    Parameters
    ----------
    requirement_strs : Tuple[str]
        Tuple of requirements strings to retrieve metadata for.

    Returns
    ----------
    OrderedDict
        Ordered dictionary lexicographically presorted on keys in ascending
        order (as detailed by the :func:`get_requirements_dict_metadata`
        function).
    '''

    # Lexicographically sorted tuple of these strings.
    requirement_strs_sorted = iterables.sort_ascending(requirement_strs)

    # List of all requirement objects parsed from these requirement strings.
    requirements = get_requirements(*requirement_strs_sorted)

    # Ordered dictionary synopsizing these requirements.
    return OrderedDict(
        (requirement.project_name, get_requirement_metadata(requirement))
        for requirement in requirements)

# ....................{ CONVERTERS ~ tuple-to-dict         }....................
@type_check
def convert_requirements_tuple_to_dict(
    requirements_tuple: SequenceTypes) -> dict:
    '''
    Convert the passed tuple of :mod:`setuptools`-specific requirements strings
    into a dictionary of such strings.

    This tuple is assumed to contain strings of the form
    ``{project_name} {project_specs}``, where:

    * ``{project_name}`` is the :mod:`setuptools`-specific project name of a
      third-party dependency (e.g., ``NetworkX``).
    * ``{project_specs}`` is a :mod:`setuptools`-specific requirements string
      constraining this dependency (e.g., ``>= 1.11``).

    Each key-value pair of the resulting dictionary maps from each such
    ``{project_name}`` to the corresponding ``{project_specs}``.

    Parameters
    ----------
    requirements_tuple : SequenceTypes
        Tuple of :mod:`setuptools`-specific requirements strings in the format
        described above.

    Returns
    ----------
    dict
        Dictionary of :mod:`setuptools`-specific requirements strings in the
        format described above.
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
    Convert the passed dictionary of :mod:`setuptools`-specific requirements
    strings into a tuple of such strings.

    This dictionary is assumed to map from the :mod:`setuptools`-specific
    project name of a third-party dependency (e.g., ``NetworkX``) to the suffix
    of a :mod:`setuptools`-specific requirements string constraining this
    dependency (e.g., ``>= 1.11``). Each element of the resulting tuple is a
    string of the form `{key} {value}`, converted from a key-value pair of this
    dictionary in arbitrary order.

    Parameters
    ----------
    requirements_dict : MappingType
        Dictionary of :mod:`setuptools`-specific requirements strings in the
        format described above.

    Returns
    ----------
    tuple
        Tuple of :mod:`setuptools`-specific requirements strings in the format
        described above.
    '''

    return convert_requirements_dict_keys_to_tuple(
        requirements_dict, *requirements_dict.keys())


@type_check
def convert_requirements_dict_keys_to_tuple(
    requirements_dict: MappingType, *requirement_names: str) -> tuple:
    '''
    Convert all key-value pairs of the passed dictionary of :mod:`setuptools`-
    specific requirements strings whose keys are the passed strings into a
    tuple of :mod:`setuptools`-specific requirements strings.

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
        Tuple of :mod:`setuptools`-specific requirements strings in the above format.

    Raises
    ----------
    BetseLibException
        If the passed key is *not* a key of this dictionary.

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
    Convert the key-value pair of the passed dictionary of :mod:`setuptools`-
    specific requirements strings whose key is the passed string into a
    :mod:`setuptools`-specific requirements string.

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
    BetseLibException
        If the passed key is *not* a key of this dictionary.

    See Also
    ----------
    :func:`convert_requirements_dict_to_tuple`
        Further details on the format of this dictionary and resulting string.
    '''

    # If this name is unrecognized, raise an exception.
    if requirement_name not in requirements_dict:
        raise BetseLibException(
            'Dependency "{}" unrecognized.'.format(requirement_name))

    # String constraining this requirement (e.g., ">= 3.14159265359", "[svg]")
    # if any or "None" otherwise.
    requirement_constraints = requirements_dict[requirement_name]

    # If this requirement is constrained, return the concatenation of this
    # requirement's name and constraints.
    if requirement_constraints:
        return '{} {}'.format(requirement_name, requirement_constraints)
    # Else, return only this requirement's name.
    else:
        return requirement_name

# ....................{ IMPORTERS                          }....................
@type_check
def import_requirements_dict_keys(
    requirements_dict: MappingType, *requirement_names: str) -> object:
    '''
    Import and return the top-level module object satisfying each
    :mod:`setuptools`-formatted requirements string produced by concatenating
    each passed key and corresponding value of the passed dictionary (e.g.,
    converting key ``six`` and value ``>= 7.9`` into the requirements string
    ``six >= 7.9``).

    Parameters
    ----------
    requirements_dict : MappingType
        Dictionary of requirements strings to import and yield a subset of.
    requirement_names : tuple[str]
        Tuple of keys identifying the key-value pairs of this dictionary to be
        imported and yielded.

    See Also
    ----------
    :func:`import_requirement_str`
        Further details.
    '''

    # Tuple of requirements strings converted by concatenating each such key and
    # value of this dictionary.
    requirement_strs = convert_requirements_dict_keys_to_tuple(
        requirements_dict, *requirement_names)

    # Validate these requirements strings.
    return import_requirement_str(*requirement_strs)


@type_check
def import_requirement_str(*requirement_strs: str) -> object:
    '''
    Import and return the top-level module object satisfying each passed
    :mod:`setuptools`-formatted requirements string (e.g., ``six >= 7.9``).

    Parameters
    ----------
    requirement_strs : tuple[str]
        Tuple of all such requirements strings to import and return modules of.

    Returns
    ----------
    object
        Either:
        * If more than one such module was imported, the tuple of these modules.
        * Else, the single imported module.

    See Also
    ----------
    :func:`import_requirement`
        Further details.
    '''

    # List of all requirement objects parsed from these requirement strings.
    requirements = get_requirements(*requirement_strs)

    # Tuple of all modules imported from these requirement objects.
    requirements_module = tuple(
        import_requirement(requirement) for requirement in requirements)

    # Return:
    #
    # * If more than one such module was imported, this tuple of these modules.
    # * Else, the single module imported above.
    return (
        requirements_module if len(requirements_module) > 1 else
        requirements_module[0])


@type_check
def import_requirement(requirement: Requirement) -> ModuleType:
    '''
    Import and return the top-level module object satisfying the passed
    :mod:`setuptools`-specific requirement.

    Parameters
    ----------
    requirement : Requirement
        Object describing this module or package's required name and version.

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
    package_name = DISTUTILS_PROJECT_NAME_TO_MODULE_NAME[
        requirement.project_name]

    # Log this importation, which can often have unexpected side effects.
    logs.log_debug('Importing third-party package "%s"...', package_name)

    # Import and return this package.
    return modules.import_module(package_name)
