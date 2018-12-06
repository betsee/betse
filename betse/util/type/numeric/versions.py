#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level versioning facilities.

Version Specifiers
----------
For generality, most callables defined by this submodule accept `PEP
440-compliant`_ version specifiers in any of several permissible formats. As
encapsulated by the :data:`VersionTypes` tuple, these are:

* :class:`str`, specifying versions in ``.``-delimited positive integer format
  (e.g., ``2.4.14.2.1.356.23``).
* :class:`tuple`, specifying versions as one or more positive integers (e.g.,
  ``(2, 4, 14, 2, 1, 356, 23)``),
* :class:`VersionSetuptoolsType`, specifying versions as instance variables
  convertible into both of the prior formats (e.g.,
  ``VersionSetuptoolsType('2.4.14.2.1.356.23')``).

.. _PEP 440-compliant:
   https://www.python.org/dev/peps/pep-0440
'''

# ....................{ IMPORTS                           }....................
import pkg_resources
from betse.util.type.types import (
    type_check, VersionSetuptoolsType, VersionTypes)

# ....................{ TESTERS ~ greater                 }....................
@type_check
def is_greater_than(version_1: VersionTypes, version_2: VersionTypes) -> bool:
    '''
    ``True`` only if the first passed version is strictly greater than the
    second passed version.

    Parameters
    ----------
    version_1 : VersionTypes
        First version (e.g., ``2.0.1``, ``(2, 0, 1)``) to be tested.
    version_2 : VersionTypes
        Second version (e.g., ``1.0.2``, ``(1, 0, 2)``) to be tested.

    Returns
    ----------
    bool
        ``True`` only if the first version is greater than the second version.

    Raises
    ----------
    pkg_resources.packaging.version.InvalidVersion
        If either version is *not* `PEP 440-compliant`_.

    .. _PEP 440-compliant:
       https://www.python.org/dev/peps/pep-0440
    '''

    # Nothing is hard. Everything is easy
    return to_comparable(version_1) > to_comparable(version_2)


@type_check
def is_greater_than_or_equal_to(
    version_1: VersionTypes, version_2: VersionTypes) -> bool:
    '''
    ``True`` only if the first passed version is greater than or equal to the
    second passed version.

    Parameters
    ----------
    version_1 : VersionTypes
        First version (e.g., ``2.0.1``, ``(2, 0, 1)``) to be tested.
    version_2 : VersionTypes
        Second version (e.g., ``1.0.2``, ``(1, 0, 2)``) to be tested.

    Returns
    ----------
    bool
        ``True`` only if the first version is greater than or equal to the
        second version.

    Raises
    ----------
    pkg_resources.packaging.version.InvalidVersion
        If either version is *not* `PEP 440-compliant`_.

    .. _PEP 440-compliant:
       https://www.python.org/dev/peps/pep-0440
    '''

    # Nothing is easy, until it is.
    return to_comparable(version_1) >= to_comparable(version_2)

# ....................{ TESTERS ~ less                    }....................
@type_check
def is_less_than(version_1: VersionTypes, version_2: VersionTypes) -> bool:
    '''
    ``True`` only if the first passed version is strictly less than the second
    passed version.

    Parameters
    ----------
    version_1 : VersionTypes
        First version (e.g., ``2.0.1``, ``(2, 0, 1)``) to be tested.
    version_2 : VersionTypes
        Second version (e.g., ``1.0.2``, ``(1, 0, 2)``) to be tested.

    Returns
    ----------
    bool
        ``True`` only if the second version is less than the first version.

    Raises
    ----------
    pkg_resources.packaging.version.InvalidVersion
        If either version is *not* `PEP 440-compliant`_.

    .. _PEP 440-compliant:
       https://www.python.org/dev/peps/pep-0440
    '''

    # Some things are hard. This is not those things.
    return to_comparable(version_1) < to_comparable(version_2)


@type_check
def is_less_than_or_equal_to(
    version_1: VersionTypes, version_2: VersionTypes) -> bool:
    '''
    ``True`` only if the first passed version is less than or equal to the
    second passed version.

    Parameters
    ----------
    version_1 : VersionTypes
        First version (e.g., ``2.0.1``, ``(2, 0, 1)``) to be tested.
    version_2 : VersionTypes
        Second version (e.g., ``1.0.2``, ``(1, 0, 2)``) to be tested.

    Returns
    ----------
    bool
        ``True`` only if the second version is less than or equal to the first
        version.

    Raises
    ----------
    pkg_resources.packaging.version.InvalidVersion
        If either version is *not* `PEP 440-compliant`_.

    .. _PEP 440-compliant:
       https://www.python.org/dev/peps/pep-0440
    '''

    # Everything is easy, until it isn't.
    return to_comparable(version_1) <= to_comparable(version_2)

# ....................{ CONVERTERS                        }....................
@type_check
def to_comparable(version: VersionTypes) -> VersionSetuptoolsType:
    '''
    Passed version converted into a **comparable version type** (i.e., type
    suitable for use both as parameters to callables accepting arbitrary
    version specifiers *and* as operands to numeric operators comparing such
    specifiers).

    Specifically, if the passed version is of type:

    * :class:`str` *or* :class:`tuple`, a new :class:`VersionSetuptoolsType`
      instance is instantiated and returned from this version.
    * :class:`VersionSetuptoolsType`, the same
      :class:`VersionSetuptoolsType` is returned as is.

    Caveats
    ----------
    The version returned by this function is *only* safely comparable with
    versions of the same type. In particular, the
    :class:`VersionSetuptoolsType` type does *not* necessarily support direct
    comparison with either the :class:`tuple` *or* `class:`str` version types;
    :class:`VersionSetuptoolsType` supported both under older but *not* newer
    versions of :mod:`setuptools`. *shakes fist*

    Parameters
    ----------
    version : VersionTypes
        Version (e.g., ``2.0.1``, ``(2, 0, 1)``) to be converted.

    Returns
    ----------
    VersionSetuptoolsType
        :mod:`setuptools`-specific version converted from this version.

    Raises
    ----------
    pkg_resources.packaging.version.InvalidVersion
        If this version is *not* `PEP 440-compliant`_.

    .. _PEP 440-compliant:
       https://www.python.org/dev/peps/pep-0440

    Examples
    ----------
    >>> from betse.util.type.numeric import versions
    >>> versions.to_comparable('2.4.14.2.1.356.23')
    Version('2.4.14.2.1.356.23')
    >>> versions.to_comparable((2, 4, 14, 2, 1, 356, 23))
    (2, 4, 14, 2, 1, 356, 23)
    '''

    # Avoid circular import dependencies.
    from betse.util.type.text import strs

    # Note that there are *MANY* approaches for converting strings into
    # comparable versions, including:
    #
    # * The packaging.version.parse() function, whose API resembles that of
    #   setuptools called below. In theory, the "packaging" package *SHOULD* be
    #   a setuptools dependency and hence always importable wherever the
    #   "pkg_resources" and "setuptools" packages are importable. In practice,
    #   this does *NOT* appear to be the case.
    # * The "distutils.version.LooseVersion" class, whose API resembles that of
    #   setuptools called below. However, whereas the latter API is PEP
    #   440-conformant, the latter API is not.
    #
    # Note also that the pkg_resources.parse_version() called below to
    # implement this conversion *ONLY*  accepts string parameters under older
    # (and possibly current) setuptools versions. Passing version tuples, in
    # particular, raises an exception resembling:
    #
    #     >>> pkg_resources.parse_version((1, 2, 3))
    #     TypeError: expected string or bytes-like object
    #
    # Setuptools. It's actually good for something. Sort of. Who knew?
    #
    # If this is a tuple-formatted version, convert this tuple into a
    # string-formatted version first. Why? Because older (and possibly current)
    # implementations of the pkg_resources.parse_version() explicitly raise
    # exceptions when passed parameters of type other than "str". *MEGASIGH*
    if isinstance(version, tuple):
        # Tuple-formatted version such that each item is guaranteed to be a
        # string rather than an integer, as required for joining these items
        # together into a composite string-formatted version.
        version_str_tuple = tuple(
            str(version_part) for version_part in version)

        # String-formatted version converted from this tuple.
        version = strs.join_on_dot(version_str_tuple)

    # This version *MUST* now be either a string or a setuptools object.
    return (
        # If this is a string-formatted version, convert this string into a
        # setuptools-specific version.
        pkg_resources.parse_version(version) if isinstance(version, str)
        # Else, this *MUST* by definition be a setuptools-specific version. In
        # this case, return this version as is.
        else version
    )
