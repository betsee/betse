#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level external process facilities.
'''

# ....................{ IMPORTS                            }....................
from betse.util.io import logs
from betse.util.path import paths
import sys

# ....................{ CONSTANTS                          }....................
EXIT_STATUS_SUCCESS = 0
'''
Exit status returned on process success.
'''

EXIT_STATUS_FAILURE_DEFAULT = 1
'''
Canonical exit status returned on process failure.
'''

# ....................{ GETTERS                            }....................
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

def exit(exit_status: int = 1, exit_message: str = '') -> None:
    '''
    Exit from the current process with the passed exit status (defaulting to 1),
    logging the passed message (defaulting to the empty string) if nonempty.

    If such exit status is 0 signifying success, such message will be logged
    with level `INFO`; else, such such message will be logged with level
    `ERROR`.
    '''
    assert isinstance(exit_message, str),\
        '"{}" not a string.'.format(exit_message)

    # Log such message if nonempty.
    if exit_message:
        # Get an appropriate child logger.
        logger = logs.get()

        # If such status signifies success, log accordingly.
        if exit_status == 0:
            logger.info(exit_message)
        # Else, such status signifies failure. Log accordingly.
        else:
            logger.error(exit_message)

    # Exit the current process with such status.
    sys.exit(exit_status)

# --------------------( WASTELANDS                         )--------------------
    # printing *without* logging the passed message (defaulting to the empty
    # string) to standard error if nonempty.
