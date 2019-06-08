#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Linux-specific facilities.
'''

# ....................{ IMPORTS                           }....................
import platform
# from betse.util.io.log import logs
from betse.util.type.decorator.decmemo import func_cached
from betse.util.type.types import StrOrNoneTypes

# ....................{ TESTERS                           }....................
@func_cached
def is_linux() -> bool:
    '''
    ``True`` only if the current platform is either Linux or a platform
    masquerading to a reasonably accurate degree as Linux (e.g., the Windows
    Subsystem for Linux (WSL) but *not* Cygwin Windows, which is sufficiently
    different to warrant differentiation).
    '''

    return platform.system() == 'Linux'


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

# ....................{ GETTERS ~ distro                  }....................
@func_cached
def get_distro_name() -> str:
    '''
    Human-readable name of the current Linux distribution.

    Specifically, this function returns:

    * If the :func:`platform.linux_distribution` function is available, the
      first element of the 3-tuple returned by this function (e.g.,
      ``CentOS``).
    * Else, the string returned by the :func:`platform.system` function.
      Since this is probably the generic string ``Linux``, this is only a
      fallback.
    '''

    # Avoid circular import dependencies.
    from betse.lib import libs
    from betse.util.io.error.errwarning import ignoring_deprecations
    from betse.util.os.brand import windows

    # Linux distribution name to be returned.
    distro_name = None

    # If the deprecated Linux-specific platform.linux_distribution()
    # function is still defined...
    if hasattr(platform, 'linux_distribution'):
        # Ignore all deprecations emitted by Python 3.5 through 3.7 concerning
        # the pending removal of the platform.linux_distribution() function...
        with ignoring_deprecations():
            # Defer to the second item of the 3-tuple "(distname, version, id)"
            # Defer to the first item of the 3-tuple "(distname, version, id)"
            # returned by this function (e.g., "CentOS"). Since
            # platform.system() returns the low-level kernel name "Linux", the
            # latter name is ignored if feasible.
            distro_name = platform.linux_distribution()[0]
    # Else if the optional third-party Linux-specific "distro" package is
    # importable, defer to this package.
    elif libs.is_runtime_optional('distro'):
        # Import this package.
        distro = libs.import_runtime_optional('distro')

        # Defer to the name reported by this package.
        distro_name = distro.name()
    # Else if this is actually the Windows Subsystem for Linux (WSL)
    # masquerading as Linux, return "Windows (WSL)" rather than "Linux".
    elif windows.is_windows_wsl():
        distro_name = 'Windows (WSL)'
    # Else, this is Python >= 3.8 *AND* the "distro" package is unimportable.
    # In this case, default to the platform name returned by the
    # platform.system() function. Since this typically corresponds to the
    # low-level name of the current kernel (e.g., "Linux") rather than the
    # high-level name of the current Linux distribution (e.g., "Arch",
    # "Gentoo"), this is only a fallback of last resort.
    else:
        distro_name = platform.system()

    # Return this version specifier.
    return distro_name


@func_cached
def get_distro_version() -> str:
    '''
    Human-readable ``.``-delimited version specifier string of the current
    platform.

    Specifically, This function returns:

    * Under Linux, the version reported by this Linux distribution as follows:

      * If the deprecated :func:`platform.linux_distribution` function is still
        defined, the second element of the 3-tuple returned by this function
        (e.g., ``6.4``).
      * Else if the third-party Linux-specific "distro" package (referenced
    # by official Python documentation as a stand-in replacement for the
    # deprecated Linux-specific platform.linux_distribution() function) is
    # importable, defer to this package.
      * Else, the string returned by the :func:`platform.release` function.
        Since this is probably the non-human-readable ``.``- and
        ``-``-delimited version specifier string of the current Linux kernel
        (e.g., ``4.1.15-gentoo-r1``), this is only a fallback.
    '''

    # Avoid circular import dependencies.
    from betse.lib import libs
    from betse.util.io.error.errwarning import ignoring_deprecations

    # Version specifier to be returned.
    distro_version = None

    # If the deprecated Linux-specific platform.linux_distribution()
    # function is still defined...
    if hasattr(platform, 'linux_distribution'):
        # Ignore all deprecations emitted by Python 3.5 through 3.7 concerning
        # the pending removal of the platform.linux_distribution() function...
        with ignoring_deprecations():
            # Defer to the second item of the 3-tuple "(distname, version, id)"
            # returned by this function (e.g., "6.4"). Since platform.release()
            # returns the low-level kernel version (e.g., "4.1.15"), the latter
            # version is ignored if feasible.
            distro_version = platform.linux_distribution()[1]
    # Else if the optional third-party Linux-specific "distro" package is
    # importable, defer to this package.
    elif libs.is_runtime_optional('distro'):
        # Import this package.
        distro = libs.import_runtime_optional('distro')

        # For disambiguity, prefer the "best" (i.e., most specific and
        # accurate) version specifier published by this Linux distribution.
        distro_version = distro.version(best=True)
    # Else, this is Python >= 3.8 *AND* the "distro" package is unimportable.
    # In this case, default to the version specifier returned by the
    # platform.release() function. Since this typically corresponds to the
    # low-level version of the current kernel (e.g., "4.1.15") rather than the
    # high-level version of the current Linux distribution (e.g., "6.4"), this
    # is a fallback of last resort.
    else:
        distro_version = platform.release()

    # Return this version specifier.
    return distro_version

# ....................{ PRIVATE ~ getters                 }....................
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
    from betse.util.os.shell import shellenv

    # If the current platform is *NOT* Linux, return "None".
    if not is_linux():
        return None
    # Else, the current platform is Linux.

    # Return the string value of the ${XDG_SESSION_TYPE} environment variable
    # if defined *OR* "None" otherwise.
    return shellenv.get_var_or_none('XDG_SESSION_TYPE')
