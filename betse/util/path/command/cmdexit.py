#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level external process facilities.
'''

# ....................{ IMPORTS                            }....................
import sys
from betse.util.io.log import logs
from betse.util.type.types import type_check

# ....................{ CONSTANTS                          }....................
SUCCESS = 0
'''
Exit status signifying a process to have terminated successfully.
'''


FAILURE_DEFAULT = 1
'''
Exit status typically signifying a process to have terminated prematurely with a
fatal error.

While any exit status in the range ``[1, 255]`` signifies failure, this exit
status is the most common and hence preferred default.
'''

# ....................{ EXCEPTIONS                         }....................
def raise_success() -> None:
    '''
    Raise the :class:`SystemExit` exception with the success exit status (i.e.,
    :data:`SUCCESS`), thus halting the current process with this status if
    uncaught.
    '''

    raise SystemExit(SUCCESS)


# ....................{ TESTERS                            }....................
@type_check
def is_success(exit_status: int) -> bool:
    '''
    ``True`` only if the passed exit status signifies success.
    '''

    return exit_status == SUCCESS


@type_check
def is_failure(exit_status: int) -> bool:
    '''
    ``True`` only if the passed exit status signifies failure.
    '''

    return exit_status != SUCCESS

# ....................{ EXITERS                            }....................
def exit_with_success(error_message: str = '') -> None:
    '''
    Halt the current process with the success exit status (i.e.,
    :data:`SUCCESS`), logging the passed error message if nonempty *or* exiting
    silently otherwise.
    '''

    exit(exit_status=FAILURE_DEFAULT, exit_message=error_message)


def exit_with_failure(error_message: str = '') -> None:
    '''
    Halt the current process with the default exit status for failure (i.e.,
    :data:`FAILURE_DEFAULT`), logging the passed error message if nonempty *or*
    exiting silently otherwise.
    '''

    exit(exit_status=FAILURE_DEFAULT, exit_message=error_message)


@type_check
def exit(exit_status: int = FAILURE_DEFAULT, exit_message: str = '') -> None:
    '''
    Halt the current process with the passed exit status, logging the passed
    message if nonempty *or* exiting silently otherwise.

    This message will be logged with level:

    * :attr:`LogLevel.INFO` if this exit status signifies success.
    * :attr:`LogLevel.ERROR` otherwise.
    '''

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
