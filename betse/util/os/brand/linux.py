#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Linux-specific facilities.
'''

# ....................{ IMPORTS                           }....................
# from betse.util.io.log import logs
from betse.util.type.decorator.decmemo import func_cached
from betse.util.type.types import StrOrNoneTypes

# ....................{ TESTERS                           }....................
@func_cached
def is_mir() -> bool:
    '''
    ``True`` only if the active Python interpreter is running under the
    Linux-specific Mir compositor, implying this process to be headfull and
    hence support both CLIs and GUIs.

    Caveats
    ----------
    For sanity, this function incorrectly assumes *all* Mir compositors to
    comply with the X Desktop Group (XDG) standard by exporting the
    ``${XDG_SESSION_TYPE}`` environment variable with a value of ``mir``. As
    this is *not* necessarily the case, this function may return false
    negatives for edge-case Mir compositors.

    See Also
    ----------
    :func:`_get_xdg_session_type_or_none`
        Further details.
    '''

    return _get_xdg_session_type_or_none() == 'mir'


@func_cached
def is_wayland() -> bool:
    '''
    ``True`` only if the active Python interpreter is running under a
    Linux-specific Wayland compositor, implying this process to be headfull and
    hence support both CLIs and GUIs.

    Caveats
    ----------
    For sanity, this function incorrectly assumes *all* Wayland compositors to
    comply with the X Desktop Group (XDG) standard by exporting the
    ``${XDG_SESSION_TYPE}`` environment variable with a value of ``wayland``.
    As this is *not* necessarily the case, this function may return false
    negatives for edge-case Wayland compositors.

    See Also
    ----------
    :func:`_get_xdg_session_type_or_none`
        Further details.
    '''

    return _get_xdg_session_type_or_none() == 'wayland'

# ....................{ GETTERS                           }....................
@func_cached
def _get_xdg_session_type_or_none() -> StrOrNoneTypes:
    '''
    Lowercase machine-readable string standardized by the X Desktop Group (XDG)
    identifying the general category of Linux-specific display manager under
    which the active Python interpreter is running if any *or* ``None``
    otherwise (e.g., if this interpreter is running headless).

    Returns
    ----------
    StrOrNoneTypes
        Either:

        * If either the current platform is not Linux *or* the current shell
          environment does not define the ``${XDG_SESSION_TYPE}`` variable
          queried by this function, ``None``. In the latter case, this commonly
          but *not* necessarily implies this interpreter to be running
          headless.
        * Else, the value of the ``${XDG_SESSION_TYPE}`` variable, guaranteed
          to be one of the following strings:

            * If this interpreter is running headless under an interactive
              pseudoterminal, ``tty``.
            * If this interpreter is running headfull under:

              * The Mir compositor, ``mir``.
              * A Wayland compositor, ``wayland``.
              * The X Window System, ``x11``.

            * If something has gone non-deterministic, ``unspecified``.

    See Also
    ----------
    https://www.freedesktop.org/software/systemd/man/sd_session_get_type.html
        Low-level PAM-specific ``sd_session_get_type(3)`` function retrieving
        the same string.
    '''

    # Avoid circular import dependencies.
    from betse.util.os import oses
    from betse.util.os.shell import shellenv

    # If the current platform is *NOT* Linux, return "None".
    if not oses.is_linux():
        return None
    # Else, the current platform is Linux.

    # Return the string value of the ${XDG_SESSION_TYPE} environment variable
    # if defined *OR* "None" otherwise.
    return shellenv.get_var_or_none('XDG_SESSION_TYPE')
