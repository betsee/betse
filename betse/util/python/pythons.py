#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
High-level Python facilities.

This module provides functionality pertaining to the active Python interpreter
as a whole.

Caveats
----------
Word size-specific functions (e.g., `is_wordsize_64()`) are generally considered
poor form. Consider calling such functions _only_ where necessary.
'''

# ....................{ IMPORTS                            }....................
from betse import metadata
from betse.util.io import loggers
from collections import OrderedDict
import platform, sys

# ....................{ INITIALIZERS                       }....................
def init() -> None:
    '''
    Validate the active Python interpreter.

    This function does _not_ validate such interpreter's version, as the
    top-level module `betse.metadata` already does so at the beginning of
    application startup. Rather, this function (in order):

    . Logs a non-fatal warning if such interpreter is _not_ 64-bit.
    '''
    if is_wordsize_32():
        loggers.log_warning(
            '32-bit Python interpreter detected. '
            '{name} will be confined to low-precision datatypes and '
            'at most 4GB of available RAM, '
            'impeding the reliability and scalability of modelling. '
            'Consider running {name} only under '
            '64-bit Python interpreters.'.format(name=metadata.NAME)
        )

# ....................{ TESTERS ~ arch                     }....................
def is_wordsize_32():
    '''
    `True` if the active Python interpreter is **32-bit** (i.e., was compiled
    with a 32-bit toolchain into a 32-bit executable).
    '''
    return not is_wordsize_64()

def is_wordsize_64():
    '''
    `True` if the active Python interpreter is **64-bit** (i.e., was compiled
    with a 64-bit toolchain into a 64-bit executable).
    '''
    # Avoid circular import dependencies.
    from betse.util.type import ints

    # There exist several alternative means of testing the same condition: e.g.,
    #
    #     return 'PROCESSOR_ARCHITEW6432' in os.environ
    #
    # The current approach, however, is more portable and hence ideal.
    return sys.maxsize > ints.INT_VALUE_MAX_32_BIT

# ....................{ GETTERS                            }....................
def get_metadata() -> OrderedDict:
    '''
    Get an ordered dictionary synopsizing the active Python interpreter.

    This function aggregates the metadata reported by the reasonably
    cross-platform module `platform` into a simple dictionary.
    '''
    # Such dictionary.
    metadata = OrderedDict((
        ('implementation', platform.python_implementation()),
        ('version', platform.python_version()),
        ('vcs revision', platform.python_revision()),
        ('vcs branch', platform.python_branch()),
        ('compiler', platform.python_compiler()),
    ))

    # 2-tuple providing such interpreter's build number and date as strings.
    python_build = platform.python_build()

    # Append such metadata.
    metadata['build number'] = python_build[0]
    metadata['build data'] = python_build[1]

    # Get such dictionary.
    return metadata

# --------------------( WASTELANDS                         )--------------------
        # Avoid circular import dependencies.
    #FUXME: Not terribly human-readable and hence hardly ideal.
    # ordered_dict[' platform.platform(aliased = True)

# def format() -> str:
#     '''
#     Get a human-readable string synopsizing the current system.
#     '''
