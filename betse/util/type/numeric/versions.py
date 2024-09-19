#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2025 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level versioning facilities.

Design
------
For generality, most callables defined by this submodule accept `PEP
440-compliant`_ version specifiers in any of several permissible formats. As
encapsulated by the :data:`VersionTypes` tuple, these are:

* :class:`str`, specifying versions in ``.``-delimited positive integer format
  (e.g., ``2.4.14.2.1.356.23``).
* :class:`tuple`, specifying versions as one or more positive integers (e.g.,
  ``(2, 4, 14, 2, 1, 356, 23)``),
* :class:`VersionSetuptoolsTypes`, specifying versions as instance variables
  convertible into both of the prior formats (e.g.,
  ``VersionSetuptoolsTypes('2.4.14.2.1.356.23')``).

.. _PEP 440-compliant:
   https://www.python.org/dev/peps/pep-0440
'''

# ....................{ IMPORTS                            }....................
from beartype.typing import Tuple
from betse.exceptions import BetseVersionException
from betse.util.type.decorator.decmemo import func_cached
from betse.util.type.types import RegexCompiledType, VersionTypes
from collections.abc import Collection as CollectionABC

# ....................{ TESTERS ~ greater                  }....................
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
    -------
    bool
        ``True`` only if the first version is greater than the second version.

    Raises
    ------
    BetseVersionException
        If either version is *not* `PEP 440-compliant`_.

    .. _PEP 440-compliant:
       https://www.python.org/dev/peps/pep-0440
    '''

    # Nothing is hard. Everything is easy.
    return to_comparable(version_1) > to_comparable(version_2)


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
    -------
    bool
        ``True`` only if the first version is greater than or equal to the
        second version.

    Raises
    ------
    BetseVersionException
        If either version is *not* `PEP 440-compliant`_.

    .. _PEP 440-compliant:
       https://www.python.org/dev/peps/pep-0440
    '''

    # Nothing is easy, until it is.
    return to_comparable(version_1) >= to_comparable(version_2)

# ....................{ TESTERS ~ less                     }....................
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
    -------
    bool
        ``True`` only if the second version is less than the first version.

    Raises
    ------
    BetseVersionException
        If either version is *not* `PEP 440-compliant`_.

    .. _PEP 440-compliant:
       https://www.python.org/dev/peps/pep-0440
    '''

    # Some things are hard. This is not those things.
    return to_comparable(version_1) < to_comparable(version_2)


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
    -------
    bool
        ``True`` only if the second version is less than or equal to the first
        version.

    Raises
    ------
    BetseVersionException
        If either version is *not* `PEP 440-compliant`_.

    .. _PEP 440-compliant:
       https://www.python.org/dev/peps/pep-0440
    '''

    # Everything is easy, until it isn't.
    return to_comparable(version_1) <= to_comparable(version_2)

# ....................{ CONVERTERS                         }....................
def to_comparable(version: VersionTypes) -> Tuple[int, ...]:
    '''
    Passed version converted into a **comparable version type** (i.e., type
    suitable for use both as parameters to callables accepting arbitrary
    version specifiers *and* as operands to numeric operators comparing such
    specifiers).

    Specifically, if the passed version is of type:

    * :class:`str`:

      * This string is munged to comply with industry-standard semantics.
        Specifically:

        * All tildes (i.e., ``~`` characters) in this string are globally
          replaced with hyphens (i.e., ``-`` characters). For unknown reasons,
          older (and possibly current) implementations of the low-level
          :func:`pkg_resources.parse_version` function underlying this
          higher-level function sort versions containing tilde characters as
          strictly less than versions *not* containing tilde characters (e.g.,
          sorting ``5.9.0~a1`` as less than ``5.7.0``). This is blatantly
          wrong, of course. Since version strings associated with dependencies
          often contain tilde characters, failure to address this subtle issue
          would induce subtle issues elsewhere throughout this application.

      * A new :class:`VersionSetuptoolsTypes` instance is instantiated and
        returned from this version.

    * :class:`tuple`, a new :class:`VersionSetuptoolsTypes` instance is
      instantiated and returned from this version.
    * :class:`VersionSetuptoolsTypes`, the same
      :class:`VersionSetuptoolsTypes` is returned as is.

    Caveats
    -------
    The version returned by this function is *only* safely comparable with
    versions of the same type. In particular, the
    :class:`VersionSetuptoolsTypes` type does *not* necessarily support direct
    comparison with either the :class:`tuple` *or* `class:`str` version types;
    :class:`VersionSetuptoolsTypes` supported both under older but *not* newer
    versions of :mod:`setuptools`. *shakes fist*

    Parameters
    ----------
    version : VersionTypes
        Version (e.g., ``2.0.1``, ``(2, 0, 1)``) to be converted.

    Returns
    -------
    VersionSetuptoolsTypes
        :mod:`setuptools`-specific version converted from this version.

    Raises
    ------
    BetseVersionException
        If this version is *not* `PEP 440-compliant`_.

    .. _PEP 440-compliant:
       https://www.python.org/dev/peps/pep-0440

    Examples
    --------
    .. code-block:: python

       >>> from betse.util.type.numeric import versions
       >>> versions.to_comparable('2.4.14.2.1.356.23')
       Version('2.4.14.2.1.356.23')
       >>> versions.to_comparable((2, 4, 14, 2, 1, 356, 23))
       (2, 4, 14, 2, 1, 356, 23)
    '''

    # ....................{ IMPORTS                        }....................
    # Avoid circular import dependencies.
    from betse.util.type.text.regexes import replace_substrs

    # ....................{ LOCALS                         }....................
    # Version passed to this function, preserved for subsequent embedding in
    # human-readable exception messages.
    version_old = version

    # ....................{ STRING                         }....................
    # Note that there are *MANY* approaches for converting strings into
    # comparable versions, including:
    # * The packaging.version.parse() function, whose API resembles that of
    #   setuptools called below. In theory, the "packaging" package *SHOULD* be
    #   a setuptools dependency and hence always importable wherever the
    #   "pkg_resources" and "setuptools" packages are importable. In practice,
    #   this does *NOT* appear to be the case.
    # * The "distutils.version.LooseVersion" class, whose API resembles that of
    #   setuptools called below. However, whereas the latter API is PEP
    #   440-conformant, the former API is not.
    # * The setuptools.pkg_resources.parse_version() function previously called
    #   below. As "setuptools" has since deprecated "pkg_resources", however,
    #   deferring to that function is no longer a reasonable proposition.

    # If this is a string-formatted version...
    if isinstance(version, str):
        # Compiled regular expression matching all invalid delimiters (e.g.,
        # "a", "b", "rc") in "."-delimited version specifiers.
        VERSION_INVALID_DELIMITER_REGEX = _get_version_invalid_delimiter_regex()

        # Sanitized version, replacing each invalid delimiter in this
        # unsanitized version with the valid delimiter ".".
        version_sanitized = replace_substrs(
            text=version,
            regex=VERSION_INVALID_DELIMITER_REGEX,
            replacement='.',
        )

        # List of one or more string-formatted version components split from
        # this sanitized version on "." delimiters. Note that this behaves as
        # expected when this version contains *NO* "." delimiters.
        version = version_sanitized.split('.')
    assert isinstance(version, CollectionABC), (
        f'{repr(version)} not version collection.')

    # ....................{ TUPLE                          }....................
    # Attempt to coerce this collection of possibly non-integer objects (e.g.,
    # string-formatted version components) into a tuple of one or more
    # integer-formatted version components.
    try:
        version = tuple(int(version_item) for version_item in version)
    except Exception as exception:
        raise BetseVersionException(
            f'PEP 440 version {repr(version_old)} invalid '
            f'(i.e., contains one or more characters that are neither '
            f'integers nor "." delimiters).'
        ) from exception

    # Return this tuple-formatted version as is.
    return version

# ....................{ PRIVATE ~ getters                  }....................
@func_cached
def _get_version_invalid_delimiter_regex() -> RegexCompiledType:
    '''
    Compiled regular expression matching all **invalid delimiters** in
    ``"."``-delimited version specifiers.

    Invalid delimiters typically denote pre-release "sub-versions," including
    such substrings as:

    * The ``"a"`` in ``"149.0a5"``.
    * The ``"b"`` in ``"1.3.0b0"``.
    * The ``"rc"`` in ``"0.19.0rc1"``.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.text.regexes import compile_regex

    # Create, return, and cache this expression.
    return compile_regex(
        # One or more instances of an invalid delimiter, defined as...
        r'(?:'
        # Dismantled, this is:
        # * "a", a single character signifying an alpha pre-release.
        # * "b", a single character signifying a beta pre-release.
        # * "~", a PySide2-specific delimiter intended to signify a "-". Look.
        #   Don't ask. Just accept that the real world is a nasty place.
        # * "-", a non-standard delimiter in some version specifiers.
        r'[ab~-]|'
        # Release candidate pre-release suffix.
        r'rc'
        r')+'
    )
