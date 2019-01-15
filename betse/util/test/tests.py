#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level **test suite** (e.g., collection of one or more functional or unit
tests collectively exercising problematic features of this application)
facilities.
'''

# ....................{ IMPORTS                           }....................
# from betse.util.io.log import logs
from betse.util.type.types import type_check

# ....................{ GLOBALS                           }....................
_IS_TESTING = False
'''
``True`` only if the active Python interpreter is running tests (e.g., with the
:mod:`pytest` test harness).
'''

# ....................{ TESTERS                           }....................
def is_testing() -> bool:
    '''
    ``True`` only if the active Python interpreter is currently running tests
    (e.g., with the :mod:`pytest` test harness).
    '''

    return _IS_TESTING

# ....................{ SETTERS                           }....................
@type_check
def set_testing(is_testing: bool) -> None:
    '''
    Set whether or not the active Python interpreter is currently running tests
    (e.g., with the :mod:`pytest` test harness).

    Parameters
    ----------
    is_testing : bool
        ``True`` only if this interpreter is currently running tests.
    '''

    # Allow this global to be locally redefined.
    global _IS_TESTING

    # Redefine this global.
    _IS_TESTING = is_testing
