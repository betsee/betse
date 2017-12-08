#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level versioning facilities.
'''

# ....................{ IMPORTS                            }....................
from betse.util.type.types import type_check
from packaging import version

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
    https://stackoverflow.com/a/11887885/2809027
        StackOverflow overflow inspiring this implementation.
    '''

    # Setuptools. It's actually good for something.
    return version.parse(version1) >= version.parse(version2)
