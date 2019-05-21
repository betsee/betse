#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
POSIX (Portable Operating System Interface)-specific facilities.
'''

# ....................{ IMPORTS                           }....................
# from betse.util.io.log import logs
from betse.util.type.decorator.decmemo import func_cached

# ....................{ TESTERS                           }....................
@func_cached
def is_x11() -> bool:
    '''
    ``True`` only if the active Python interpreter is running under the X
    Window System (hereafter X11), implying this process to be headfull and
    hence support both CLIs and GUIs.

    Caveats
    ----------
    This function returning ``True`` does *not* necessarily imply the current
    platform to be a Linux distribution. Unlike most Linux-centric protocols,
    support for the X11 protocol is sufficiently widespread across non-Linux
    POSIX-compatible platforms as to be effectively cross-platform. Unlike the
    Linux-specific :func:`is_mir` and :func:`is_wayland` testers, this
    Linux-agnostic tester detects X11 by generically detecting the X11-specific
    ``${DISPLAY}`` environment variable rather than the Linux-specific
    ``${XDG_SESSION_TYPE}`` environment variable if any.
    '''

    # Avoid circular import dependencies.
    from betse.util.os import oses
    from betse.util.os.shell import shellenv

    # If the current platform is *NOT* POSIX-compatible, return false.
    if not oses.is_posix():
        return False
    # Else, the current platform is POSIX-compatible.

    # Return true only if the active Python interpreter inherited the
    # X11-specific ${DISPLAY} variable from its parent shell environment.
    return shellenv.is_var('DISPLAY')
