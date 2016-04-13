#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
High-level Python facilities.

This module provides functionality pertaining to the active Python interpreter
as a whole.

Caveats
----------
Word size-specific functions (e.g., `is_wordsize_64()`) are generally considered
poor form. Call these functions _only_ where necessary.
'''

# ....................{ IMPORTS                            }....................
from betse import metadata
from betse.util.io import logs
from collections import OrderedDict
import platform, sys

# ....................{ INITIALIZERS                       }....................
def init() -> None:
    '''
    Validate the active Python interpreter.

    This function does _not_ validate this interpreter's version, as the
    top-level module `betse.metadata` already does so at the beginning of
    application startup. Rather, this function (in order):

    . Logs a non-fatal warning if this interpreter is _not_ 64-bit.
    '''

    # If this Python interpreter is 32- rather than 64-bit, log a non-fatal
    # warning. While technically feasible, running BETSE under 32-bit Python
    # interpreters imposes non-trivial constraints detrimental to sanity.
    if is_wordsize_32():
        logs.log_warning(
            '32-bit Python interpreter detected. '
            '{name} will be confined to low-precision datatypes and '
            'at most 4GB of available RAM, '
            'impeding the reliability and scalability of modelling. '
            'Consider running {name} only under '
            '64-bit Python interpreters.'.format(name=metadata.NAME)
        )

# ....................{ TESTERS                            }....................
def is_frozen() -> bool:
    '''
    `True` only if the active Python interpreter is **frozen** (i.e., embedded
    in a platform-specific compressed executable archiving `betse` and all
    transitive dependencies thereof).
    '''

    # If the "sys" module has an attribute:
    #
    # * "_MEIPASS", this is a binary frozen by PyInstaller.
    # * "frozen", this is a binary frozen by a non-PyInstaller freezer (e.g.,
    #   "py2app", "py2exe").
    return hasattr(sys, '_MEIPASS') or hasattr(sys, 'frozen')

# ....................{ TESTERS ~ arch                     }....................
def is_wordsize_32() -> bool:
    '''
    `True` only if the active Python interpreter is **32-bit** (i.e., was
    compiled with a 32-bit toolchain into a 32-bit executable).
    '''
    return not is_wordsize_64()


def is_wordsize_64() -> bool:
    '''
    `True` only if the active Python interpreter is **64-bit** (i.e., was
    compiled with a 64-bit toolchain into a 64-bit executable).
    '''

    # Avoid circular import dependencies.
    from betse.util.type import ints

    # There exist several alternative means of testing the same condition: e.g.,
    #
    #     return 'PROCESSOR_ARCHITEW6432' in os.environ
    #
    # The current approach, however, is the most portable and hence ideal.
    return sys.maxsize > ints.INT_VALUE_MAX_32_BIT

# ....................{ GETTERS                            }....................
def get_name() -> str:
    '''
    Human-readable name of the active Python interpreter's implementation (e.g.,
    `CPython`, `PyPy`).
    '''
    return platform.python_implementation()


def get_version() -> str:
    '''
    Human-readable `.`-delimited version specifier string of the active Python
    interpreter (e.g., `2.7.10`, `3.4.1`).
    '''
    return platform.python_version()

# ....................{ GETTERS ~ metadata                 }....................
def get_metadata() -> OrderedDict:
    '''
    Get an ordered dictionary synopsizing the active Python interpreter.

    This function aggregates the metadata reported by the reasonably
    cross-platform module `platform` into a simple dictionary.
    '''
    # Such dictionary.
    metadata = OrderedDict((
        ('type', get_name()),
        ('version', get_version()),
        ('vcs revision', platform.python_revision()),
        ('vcs branch', platform.python_branch()),
        ('compiler', platform.python_compiler()),
    ))

    # 2-tuple providing this interpreter's build number and date as strings.
    python_build = platform.python_build()

    # Append such metadata.
    metadata['build number'] = python_build[0]
    metadata['build data'] = python_build[1]

    # Return this dictionary.
    return metadata
