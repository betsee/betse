#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level **venv** (i.e., lightweight virtual environment containing a
distinct Python binary and site directories, optionally isolated from system
site directories) facilities.
'''

# ....................{ IMPORTS                           }....................
import sys
# from betse.util.io.log import logs
from betse.util.type.decorator.decmemo import func_cached
# from betse.util.type.types import type_check

# ....................{ TESTERS                           }....................
#FIXME: Generalize this tester to detect Anaconda-specific venvs as well. See:
#    https://stackoverflow.com/a/40099080/2809027
@func_cached
def is_venv() -> bool:
    '''
    ``True`` only if the active Python interpreter is isolated to a venv.

    Specifically, this function returns ``True`` only if this interpreter is
    isolated to a venv produced by either:

    * The official :mod:`venv` package bundled with Python >= 3.3.
    * The third-party :mod:`virtualenv` package supporting both Python 2 and 3.

    See Also
    ----------
    https://stackoverflow.com/a/42580137/2809027
        StackOverflow answer strongly inspiring this implementation.
    '''

    # Return true if this interpreter is isolated to a venv produced by
    # either...
    return (
        # "virtualenv", which uniquely defines the "sys.real_prefix"
        # attribute to the absolute dirname of the top-level directory
        # containing the system-wide Python interpreter *OR*...
        hasattr(sys, 'real_prefix') or

        # "venv", which (possibly non-uniquely) sets:
        #
        # * The "sys.base_prefix" attribute to the absolute dirname of the
        #   top-level directory containing the system-wide Python interpreter.
        # * The "sys.prefix" attribute to the absolute dirname of the
        #   top-level directory containing the venv-specific Python interpreter
        #   if any *OR* the system-wide Python interpreter otherwise.
        #
        # Note that, as Python >= 3.3 *ALWAYS* defines the "sys.base_prefix"
        # attribute, testing this attribute's existence is unnecessary.
        sys.prefix != sys.base_prefix
    )

# ....................{ GETTERS                           }....................
#FIXME: Generalize this tester to support Anaconda-specific venvs as well.
@func_cached
def get_system_prefix() -> bool:
    '''
    Absolute dirname of the top-level directory containing the system-wide
    Python interpreter, regardless of whether the active Python interpreter is
    isolated to a venv.
    '''

    # Return either...
    return (
        # The "virtualenv"-specific attribute providing this dirname if the
        # the active Python interpreter is isolated to a venv *OR*...
        sys.real_prefix
        if hasattr(sys, 'real_prefix') else
        # The standard attribute providing this dirname otherwise.
        sys.base_prefix
    )
