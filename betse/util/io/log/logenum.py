#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level logging enumerations.
'''

# ....................{ IMPORTS                            }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To avoid circular import dependencies, avoid importing from *ANY*
# application-specific modules at the top-level -- excluding those explicitly
# known *NOT* to import from this module. Since all application-specific modules
# must *ALWAYS* be able to safely import from this module at any level, these
# circularities are best avoided here rather than elsewhere.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

import enum, logging
from enum import IntEnum

# ....................{ ENUMS ~ level                      }....................
@enum.unique
class LogLevel(IntEnum):
    '''
    Enumeration of all possible **logging levels** (i.e., integer constants
    defined by the standard :mod:`logging` module, comparable according to the
    established semantics of the ``<`` comparator).

    This enumeration corresponds exactly to the ``--log-level`` CLI option.

    Comparison
    ----------
    Enumeration members are integer constants defined in increasing order.
    Enumeration members assigned smaller integers are more inclusive (i.e.,
    correspond to logging levels that log strictly more log messages than)
    enumeration members assigned larger integers: e.g.,

        # "DEBUG" is less and hence more inclusive than "INFO".
        >>> LogLevel.DEBUG < LogLevel.INFO
        True

    Usability
    ----------
    Enumeration members are integer constants and hence implicitly usable
    wherever an integer r-value (i.e., value in the right-hand side of a
    variable assignment statement) is expected: e.g.,

        # Instruct the root logger to entertain all log requests. Neat-o, eh?
        >>> import logging
        >>> logging.getLogger('').setLevel(LogLevel.ALL)
    '''


    ALL = logging.NOTSET
    '''
    Logging level signifying all messages to be logged.

    Since the :mod:`logging` module defines no constants encapsulating the
    concept of "all" (i.e., of logging everything), this is an ad-hoc constant
    expected to be smaller than the smallest constant defined by that module.
    '''


    DEBUG = logging.DEBUG
    '''
    Logging level suitable for debugging messages.
    '''


    INFO = logging.INFO
    '''
    Logging level suitable for informational messages.
    '''


    WARNING = logging.WARNING
    '''
    Logging level suitable for warning messages.
    '''


    ERROR = logging.ERROR
    '''
    Logging level suitable for error messages.
    '''


    CRITICAL = logging.CRITICAL
    '''
    Logging level suitable for critical messages.
    '''


    NONE = logging.CRITICAL + 1024
    '''
    Logging level signifying no messages to be logged.

    Since the :mod:`logging` module defines no constants encapsulating the
    concept of "none" (i.e., of logging nothing), this is an ad-hoc constant
    expected to be larger than the largest constant defined by that module.
    '''
