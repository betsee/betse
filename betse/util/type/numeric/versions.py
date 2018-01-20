#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level versioning facilities.
'''

# ....................{ IMPORTS                            }....................
import pkg_resources
from betse.util.type.types import type_check

# ....................{ TESTERS                            }....................
@type_check
def is_at_least(version1: str, version2: str) -> bool:
    '''
    ``True`` only if the version described by the first passed ``.``-delimited
    version specifier is at least as recent as the version described by the
    second such specifier.

    Parameters
    ----------
    version1: str
        First ``.``-delimited version specifier (e.g., ``2.0.1``) to test.
    version2: str
        Second ``.``-delimited version specifier (e.g., ``1.0.2``) to test.

    Returns
    ----------
    bool
        ``True`` only if the first version is as recent as the second version.

    See Also
    ----------
    https://stackoverflow.com/a/21065570/2809027
        StackOverflow overflow inspiring this implementation.
    '''

    # Setuptools. It's actually good for something.
    #
    # Note that there are *MANY* varying approaches for comparing version
    # specifiers from Python, including:
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
    # Hence, the current approach.
    return (
        pkg_resources.parse_version(version1) >=
        pkg_resources.parse_version(version2)
    )
