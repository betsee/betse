#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level external process facilities.
'''

# ....................{ IMPORTS                            }....................
import sys

from betse.util.io.log import logs
from betse.util.type import types

# ....................{ CONSTANTS                          }....................
SUCCESS = 0
'''
Exit status signifying a process to have terminated successfully.
'''

FAILURE_DEFAULT = 1
'''
Exit status typically signifying a process to have terminated prematurely with a
fatal error.

While any exit status in the range `[1, 255]` signifies failure, this exit
status is the most common and hence preferred default.
'''

# ....................{ TESTERS                            }....................
def is_success(exit_status: int) -> bool:
    '''
    `True` only if the passed exit status signifies success.
    '''
    assert types.is_int(exit_status), types.assert_not_int(exit_status)
    return exit_status == SUCCESS


def is_failure(exit_status: int) -> bool:
    '''
    `True` only if the passed exit status signifies failure.
    '''
    assert types.is_int(exit_status), types.assert_not_int(exit_status)
    return exit_status != SUCCESS

# ....................{ EXITERS                            }....................
def exit_with_failure(error_message: str = '') -> None:
    '''
    Exit from the current process with the default exit status for failure
    (i.e., `FAILURE_DEFAULT`), logging the passed error message if
    nonempty _or_ exiting silently otherwise.
    '''

    exit(exit_status=FAILURE_DEFAULT, exit_message=error_message)


def exit(
    exit_status: int = FAILURE_DEFAULT,
    exit_message: str = ''
) -> None:
    '''
    Exit from the current process with the passed exit status, logging the
    passed message if nonempty _or_ exiting silently otherwise.

    This message will be logged with level:

    * `INFO` if this exit status signifies success.
    * `ERROR` otherwise.
    '''
    assert types.is_str(exit_message), (
        types.assort_not_str(exit_message, 'Exit message'))

    # Log this message if nonempty.
    if exit_message:
        # If this status signifies success, log accordingly.
        if exit_status == SUCCESS:
            logs.log_info(exit_message)
        # Else, this status signifies failure. Log accordingly.
        else:
            logs.log_error(exit_message)

    # Exit the current process with this status.
    sys.exit(exit_status)
