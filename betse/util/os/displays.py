#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level operating system (OS)-specific display facilities.
'''

#FIXME: Call the set_headless() function if the "--headless" option is passed.

#FIXME: Submit a Stackoverflow answer encapsulating this logic. The
#osx.is_aqua() function in particular would be useful to a wide audience.

# ....................{ IMPORTS                           }....................
from betse.util.io.log import logs
from betse.util.type.decorator.decmemo import func_cached
from betse.util.type.iterable.mapping.mapcls import OrderedArgsDict
from betse.util.type.types import type_check

# ....................{ GLOBALS                           }....................
_is_headless_forced = None
'''
``True`` or ``False`` only if the :func:`set_headless` function (which
overrides the detection performed by the :func:`is_headless` function of
whether the active Python interpreter is running headless) has been called at
least once, in which case all subsequent calls to the :func:`is_headless`
function returns this boolean rather than performing that detection.

Defaults to ``None``, in which case the :func:`is_headless` function performs
such detection rather than returning the value of this boolean.
'''

# ....................{ TESTERS ~ head(full|less)         }....................
# Note that the following public testers *CANNOT* be memoized (e.g., via the
# @func_cached decorator), as the set_headless() function allows callers to
# externally modify the return values of these testers at any time. Instead,
# the private _is_headless() tester (to which these public testers defer)
# internally caches the costly detection it performs on the first call to that
# tester and then subsequently returns that cached boolean.

def is_headfull() -> bool:
    '''
    ``True`` only if the active Python interpreter is running **headfull**
    (i.e., with access to a GUI display, the common case when running under a
    conventional desktop, laptop, or tablet device).
    '''

    return not is_headless()     # Makes sense.


def is_headless() -> bool:
    '''
    ``True`` only if the active Python interpreter is running **headless**
    (i.e., with *no* access to a GUI display, often due to running remotely
    over an SSH-encrypted connection supporting only CLI input and output).

    See Also
    ----------
    :func:`_is_headless`
        Further details.
    '''

    # If the set_headless() function has been called at least once, return the
    # boolean passed to the most recent call to that function.
    if _is_headless_forced is not None:
        return _is_headless_forced
    # Else, the set_headless() function has *NOT* been called at least once. In
    # this case, intelligently detect whether this process is headless or not.

    # Return true only if this interpreter is running headless.
    return _is_headless()


@func_cached
def _is_headless() -> bool:
    '''
    ``True`` only if the active Python interpreter is running **headless**
    (i.e., with *no* access to a GUI display, often due to running remotely
    over an SSH-encrypted connection supporting only CLI input and output).

    Specifically, this function returns:

    * If the :func:`set_headless` function has been called at least once, the
      boolean passed to the most recent call to that function. Ergo, that
      function overrides the intelligent detection performed by this function.
    * Else, ``True`` only if none of the following conditions apply:

      * The current platform is Microsoft Windows. While a small subset of
        server-specific variants of Windows can and often are run headless
        (e.g., Windows Nano Server), there appears to be no known means of
        reliably distinguishing a headless from headfull Windows environment in
        pure Python. For safety, we assume the latter.
      * The current platform is not Microsoft Windows (and hence is
        POSIX-compatible) *and* either:

        * The ``${DISPLAY}`` environment variable required by the
          POSIX-compatible X11 display server is set under the current
          environment (typically due to being inherited from the parent
          environment). Note that *all* POSIX-compatible platforms of interest
          including macOS support this server.
        * The current platform is macOS *and* the :func:`macos.is_aqua` returns
          ``True``.
        * The current platform is Linux *and* either:

          * The ``${MIR_SOCKET}`` environment variable required by the
            Linux-specific Mir display server is set under this environment.
          * The ``${WAYLAND_DISPLAY}`` environment variable required by the
            Linux-specific Wayland display server is set under this
            environment.
    '''

    # Avoid circular import dependencies.
    from betse.util.os import oses
    from betse.util.os.brand import macos
    from betse.util.os.shell import shellenv

    # The active Python interpreter is headfull if and only if either...
    is_os_headfull = (
        # This is Windows, in which case this interpreter is usually headfull.
        # While certain server-specific variants of Windows can and often are
        # run headless (e.g., Windows Nano Server), there appears to be no
        # known means of reliably distinguishing a headless from headfull
        # Windows environment in pure Python. For safety, assume the latter.
        oses.is_windows() or

        # Else, this is a POSIX-compatible platform.
        #
        # Since all POSIX-compatible platforms of interest support the headfull
        # X11 display server, we efficiently test for the accessibility of this
        # server via the ${DISPLAY} environment variable inherited from the
        # parent shell environment first.
        shellenv.is_var('DISPLAY') or

        #FIXME: Unify this test with the is_linux_wayland() function, which
        #appears to be considerably more robust than the test performed here.

        # Else, all possible alternative display servers specific to the
        # current platform *MUST* be iteratively tested for.
        #
        # If this is Linux, the only remaining display servers are:
        #
        # * Mir, accessible via the ${MIR_SOCKET} environment variable.
        # * Wayland, accessible via the ${WAYLAND_DISPLAY} environment
        #   variable.
        #
        # Ergo, the current process is headfull if and only if one of these
        # variables is inherited from the parent shell environment.
        (oses.is_linux() and
         shellenv.is_var('MIR_SOCKET', 'WAYLAND_DISPLAY',)) or

        # If this is OS X, the only remaining display server is Aqua.
        (oses.is_macos() and macos.is_aqua())

        # Else, this platform is unrecognized. For safety, this platform is
        # assumed to be headless.
    )

    # Return true only if this interpreter is *NOT* headfull.
    #
    # Note that the "is_os_headfull" boolean intentionally contains the core
    # detection logic, as detecting headfull environments is fundamentally
    # more intuitive than detecting the converse.
    return not is_os_headfull

# ....................{ TESTERS ~ linux                   }....................
#FIXME: Shift into the "betse.util.os.brand.linux" submodule.
@func_cached
def is_linux_wayland() -> bool:
    '''
    ``True`` only if the active Python interpreter is running under a Wayland
    compositor-enabled Linux distribution.

    Caveats
    ----------
    For sanity, this function incorrectly assumes *all* Wayland compositors to
    comply with the X Desktop Group (XDG) standard by exporting the
    ``${XDG_SESSION_TYPE}`` environment variable with a value of ``wayland``.
    As this is *not* necessarily the case, this function may return false
    negatives for edge-case Wayland compositors.
    '''

    # Avoid circular import dependencies.
    from betse.util.os import oses
    from betse.util.os.shell import shellenv

    # If the current platform is *NOT* Linux, return False.
    if not oses.is_linux():
        return False
    # Else, the current platform is Linux.

    # String value of the ${XDG_SESSION_TYPE} environment variable if defined
    # *OR* "None" otherwise.
    xdg_session_type = shellenv.get_var_or_none('XDG_SESSION_TYPE')

    # Return True only if this value is that of a Wayland compositor.
    return xdg_session_type == 'wayland'

# ....................{ SETTERS                           }....................
@type_check
def set_headless(is_headless: bool) -> None:
    '''
    Set whether the active Python interpreter is running **headless** (i.e.,
    with *no* access to a GUI display) or not, explicitly overriding the
    implicit detection performed by the :func:`is_headfull` function of whether
    this interpreter actually is running headless or not.

    Parameters
    ----------
    is_headless : bool
        ``True`` only if the active Python interpreter is to be treated as if
        it were effectively running headless - regardless of whether it is.
    '''

    # Enable this global to be locally set.
    global _is_headless_forced

    # Log this coercion.
    logs.log_debug(
        'Coercing headless environment detection to "%r"...', is_headless)

    # Set this global to this boolean.
    _is_headless_forced = is_headless

# ....................{ GETTERS ~ metadata                }....................
def get_metadata() -> OrderedArgsDict:
    '''
    Ordered dictionary synopsizing the current display.
    '''

    # Avoid circular import dependencies.
    from betse.util.os.brand import macos

    # Return this dictionary.
    return OrderedArgsDict(
        'headless', is_headless(),
        'aqua',     macos.is_aqua(),
        'wayland',  is_linux_wayland(),
    )
