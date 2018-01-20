#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level logging filter subclasses.
'''

# ....................{ IMPORTS                            }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To avoid circular import dependencies, avoid importing from *ANY*
# application-specific modules at the top-level -- excluding those explicitly
# known *NOT* to import from this module. Since all application-specific modules
# must *ALWAYS* be able to safely import from this module at any level, these
# circularities are best avoided here rather than elsewhere.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

from betse import metadata
from betse.util.io.log.logenum import LogLevel
from betse.util.type.types import type_check
from logging import Filter, LogRecord

# ....................{ CLASSES                            }....................
class LogFilterThirdPartyDebug(Filter):
    '''
    Log filter ignoring all log records with logging levels less than or equal
    to :attr:`LogLevel.DEBUG` *and* names not prefixed by ``betse``.

    Equivalently, this log filter *only* retains log records with either:

    * Logging levels greater than :attr:`LogLevel.DEBUG`.
    * Names prefixed by ``betse``, including both:
      * ``betse``, the top-level package for BETSE.
      * ``betsee``, the top-level package for BETSEE.

    This log filter prevents ignorable debug messages logged by third-party
    frameworks (e.g., Pillow) from polluting this application's debug output.
    '''

    @type_check
    def filter(self, log_record: LogRecord) -> bool:
        '''
        ``True`` only if the passed log record is to be retained.
        '''

        # print('log record name: {}'.format(log_record.name))
        return (
            log_record.levelno > LogLevel.DEBUG or
            log_record.name.startswith(metadata.PACKAGE_NAME))


class LogFilterMoreThanInfo(Filter):
    '''
    Log filter ignoring all log records with logging levels greater than
    :attr:`LogLevel.INFO``.

    Equivalently, this log filter *only* retains log records with logging levels
    less than or equal to :attr:`LogLevel.INFO``.
    '''

    @type_check
    def filter(self, log_record: LogRecord) -> bool:
        '''
        ``True`` only if the passed log record is to be retained.
        '''

        return log_record.levelno <= LogLevel.INFO
