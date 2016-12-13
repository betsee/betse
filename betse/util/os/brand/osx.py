#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Apple OS X-specific facilities.
'''

# ....................{ IMPORTS                            }....................
from betse.exceptions import BetseOSException
from ctypes import cdll

# ....................{ CONSTANTS                          }....................
_SESSION_HAS_GRAPHIC_ACCESS = 0x0010
'''
Magic hexadecimal number defined by the OS X-specific
`/System/Library/Frameworks/Security.Framework/Headers/AuthSession.h` C header
'''

# ....................{ EXCEPTIONS                         }....................
def die_unless_os_x() -> None:
    '''
    Raise an exception unless the current platform is Apple OS X.
    '''

    # Avoid circular import dependencies.
    from betse.util.os import oses

    # If the current platform is *NOT* OS X, raise an exception.
    if not oses.is_os_x():
        raise BetseOSException(
            'Current platform {} not OS X.'.format(oses.get_name()))

# ....................{ TESTERS                            }....................
#FIXME: Implement us up.
def is_aqua() -> bool:
    '''
    `True` only if the current process is running under and hence has access to
    the OS X-specific Aqua display server.
    '''

    # Raise an exception unless the current platform is OS X.
    die_unless_os_x()
