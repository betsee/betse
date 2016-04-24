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
from betse.util.path import paths
from betse.util.type import types

# ....................{ CONSTANTS                          }....................
EXIT_STATUS_SUCCESS = 0
'''
Exit status returned on process success.
'''

EXIT_STATUS_FAILURE_DEFAULT = 1
'''
Canonical exit status returned on process failure.
'''

# ....................{ GETTERS ~ path                     }....................
def get_current_basename() -> str:
    '''
    Get the basename of the executable originating the current process.
    '''

    # Since "sys.argv[0]" is either an absolute or relative path, get only such
    # path's basename.
    return paths.get_basename(sys.argv[0])

# ....................{ EXITERS                            }....................
def exit_with_failure(exit_message: str = '') -> None:
    '''
    Exit from the current process with canonical failure exit status 1,
    logging the passed message (defaulting to the empty string) with logging
    level `ERROR` if such message is nonempty.
    '''
    exit(
        exit_status = EXIT_STATUS_FAILURE_DEFAULT,
        exit_message = exit_message,
    )


def exit(
    exit_status: int = EXIT_STATUS_FAILURE_DEFAULT,
    exit_message: str = ''
) -> None:
    '''
    Exit from the current process with the passed exit status (defaulting to
    `EXIT_STATUS_FAILURE_DEFAULT`), logging the passed message (defaulting to
    the empty string) if nonempty.

    This message will be logged with level:

    * `INFO` if this exit status signifies success.
    * `ERROR` otherwise.
    '''
    assert types.is_str(exit_message), (
        types.assort_not_str(exit_message, 'Exit message'))

    # Log this message if nonempty.
    if exit_message:
        # If this status signifies success, log accordingly.
        if exit_status == EXIT_STATUS_SUCCESS:
            logs.log_info(exit_message)
        # Else, this status signifies failure. Log accordingly.
        else:
            logs.log_error(exit_message)

    # Exit the current process with such status.
    sys.exit(exit_status)
