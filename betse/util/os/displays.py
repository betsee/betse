#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level operating system (OS)-specific display facilities.
'''

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

# ....................{ TESTERS                           }....................
@func_cached
def is_dpi_scaling() -> bool:
    '''
    ``True`` only if the current platform natively supports high-DPI scaling.

    Specifically, if the current platform is:

    * Linux *and* the active Python interpreter is running under a Wayland
      compositor, this function returns ``True``.
    * macOS, this function returns ``True``.
    * Windows >= 10, this function returns ``True``. Technically, different
      versions of Windows 10 allow end users to conditionally disable high-DPI
      scaling for different use cases. Since detecting these versions and
      settings in pure-Python is effectively infeasible, this function
      simplistically assumes all installations of Windows 10 and newer to
      unconditionally enable high-DPI scaling.

    All other platforms are assumed to *not* natively support high-DPI scaling.
    This includes both the Linux-centric X11 display server *and* all versions
    of Windows preceding Windows 10 (i.e., Windows <= 8).
    '''

    # Avoid circular import dependencies.
    from betse.util.os.brand import linux, macos, windows

    # Return true only if this platform is either...
    return (
        # A Linux distribution running Wayland *OR*...
        linux.is_wayland() or
        # macOS *OR*...
        macos.is_aqua() or
        # Windows >= 10.
        windows.is_version_10_or_newer()
    )

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
    from betse.util.os.brand import linux, macos, posix, windows

    # The active Python interpreter is headfull if and only if either...
    is_os_headfull = (
        # This is Windows, in which case this interpreter is usually headfull.
        # While certain server-specific variants of Windows can and often are
        # run headless (e.g., Windows Nano Server), there appears to be no
        # known means of reliably distinguishing a headless from headfull
        # Windows environment in pure Python. For safety, assume the latter.
        windows.is_windows() or

        # Else, this is a POSIX-compatible platform.
        #
        # Since all POSIX-compatible platforms of interest support the popular
        # X11 display server, detect this server first.
        posix.is_x11() or

        # Else, all possible alternative display servers specific to the
        # current platform *MUST* be iteratively tested for.
        #
        # If Linux, the only remaining display servers are Mir and Wayland.
        (linux.is_linux() and (linux.is_wayland() or linux.is_mir())) or

        # If macOS, the only remaining display server is Aqua.
        (macos.is_macos() and macos.is_aqua())

        # Else, this platform is unrecognized. For safety, this platform is
        # assumed to be headless.
    )

    # Return true only if this interpreter is *NOT* headfull.
    #
    # Note that the "is_os_headfull" boolean intentionally contains the core
    # detection logic, as detecting headfull environments is fundamentally
    # more intuitive than detecting the converse.
    return not is_os_headfull

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

    # If coercing headless operation, log this coercion.
    if is_headless:
        logs.log_debug('Forcing headless operation...')
    # Else, headfull operation is being coerced. In this case...
    else:
        # Log this coercion.
        logs.log_debug('Forcing headfull operation...')

        # If the current environment is detected to be headless, log a
        # non-fatal warning. While an exception could also be raised, our
        # detection heuristic is known to be imperfect.
        if _is_headless():
            logs.log_warning(
                'Headless environment detected! '
                'Forcing headfull operation under a headless environment '
                'typically raises silent segmentation faults '
                'and hence is unsupported.'
            )

    # Set this global to this boolean.
    _is_headless_forced = is_headless

# ....................{ GETTERS ~ metadata                }....................
def get_metadata() -> OrderedArgsDict:
    '''
    Ordered dictionary synopsizing the current display.
    '''

    # Avoid circular import dependencies.
    from betse.util.os.brand import linux, macos, posix

    # Return this dictionary.
    return OrderedArgsDict(
        'headless',    is_headless(),
        'dpi scaling', is_dpi_scaling(),
        'aqua',        macos.is_aqua(),
        'mir',         linux.is_mir(),
        'wayland',     linux.is_wayland(),
        'x11',         posix.is_x11(),
    )
