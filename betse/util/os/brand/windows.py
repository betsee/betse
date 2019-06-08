#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Microsoft Windows-specific facilities.

Caveats
----------
Operating system-specific logic is poor form and should be leveraged only where
necessary.
'''

# ....................{ IMPORTS                           }....................
import platform, sys
from betse.exceptions import BetseOSException
from betse.util.type.decorator.decmemo import func_cached
from betse.util.type.types import type_check, NoneType

# ....................{ TYPES                             }....................
WindowsErrorType = NoneType
'''
Windows-specific :class:`WindowsError` class if the current platform is Windows
*or* the :class:`NoneType` class otherwise.

Since all exceptions subclass the root :class:`Exception` superclass *and*
since the :class:`NoneType` class does not do so, this class is only an
exception class under Windows. In particular, attempting to either catch this
class as an exception *or* to test whether a caught exception is an instance of
this class is guaranteed to both be safe and behave as expected regardless of
the current platform. In short: don't worry, be emoji.
'''

# For use in type checking, this class is safely defined at the top-level in a
# platform-agnostic manner.
try:
    WindowsErrorType = WindowsError
except:
    pass

# ....................{ CONSTANTS ~ error codes           }....................
# For conformance, the names of all error code constants defined below are
# exactly as specified by Microsoft itself. Sadly, Python fails to provide
# these magic numbers for us.

_ERROR_INVALID_NAME = 123
'''
Microsoft Windows-specific error code indicating an invalid pathname.

See Also
----------
https://msdn.microsoft.com/en-us/library/windows/desktop/ms681382%28v=vs.85%29.aspx
    Official listing of all such codes.
'''

# ....................{ EXCEPTIONS                        }....................
def die_unless_windows() -> None:
    '''
    Raise an exception unless the current platform is Microsoft Windows.

    See Also
    ----------
    :func:`is_windows`
        Further details.
    '''

    # Avoid circular import dependencies.
    from betse.util.os import oses

    # If the current platform is *NOT* Windows, raise an exception.
    if not is_windows():
        raise BetseOSException(
            'Current platform {} not Windows.'.format(oses.get_name()))

# ....................{ TESTERS                           }....................
@func_cached
def is_windows() -> bool:
    '''
    ``True`` only if the current platform is Microsoft Windows.

    This function reports ``True`` for both vanilla and Cygwin Microsoft
    Windows (both of which commonly require special Windows-specific handling)
    but *not* the Windows Subsystem for Linux (WSL) (which accurately
    masquerades as Linux and hence does *not* commonly require special
    Windows-specific handling).
    '''

    return is_windows_vanilla() or is_windows_cygwin()


@func_cached
def is_windows_cygwin() -> bool:
    '''
    ``True`` only if the current platform is **Cygwin Microsoft
    Windows** (i.e., running the Cygwin POSIX compatibility layer).
    '''

    return sys.platform == 'cygwin'


@func_cached
def is_windows_vanilla() -> bool:
    '''
    ``True`` only if the current platform is **vanilla Microsoft
    Windows** (i.e., *not* running the Cygwin POSIX compatibility layer).
    '''

    return sys.platform == 'win32'


@func_cached
def is_windows_wsl() -> bool:
    '''
    ``True`` only if the current platform is **Windows Subsystem for
    Linux (WSL)** (i.e., the Microsoft-flavoured Linux kernel optionally
    supported by Windows 10).

    See Also
    ----------
    https://www.reddit.com/r/bashonubuntuonwindows/comments/85jghk/how_to_allow_python_to_know_im_on_windows/dvxwy9t
        Reddit post strongly inspiring this implementation.
    '''

    # Avoid circular import dependencies.
    from betse.util.os.brand import linux

    # If the active Python interpreter is *NOT* operating under a Linux kernel,
    # return false immediately.
    if not linux.is_linux():
        return False
    # Else, this interpreter is operating under a Linux kernel.

    # Flavour of this Linux kernel.
    kernel_flavour = platform.uname()[3]

    # Return true only if this is a Microsoft-flavoured Linux kernel.
    return 'microsoft' in kernel_flavour

# ....................{ TESTERS ~ path                    }....................
@type_check
def is_exception_pathname_invalid(exception: WindowsErrorType) -> bool:
    '''
    ``True`` only if the passed Windows-specific exception was raised from an
    erroneous attempt to read or write an invalid pathname.
    '''

    # If the current platform is *NOT* Windows, raise an exception.
    die_unless_windows()

    # The above type checking guarantees this exception to be an instance
    # of the Windows-specific "WindowsError" class. Test the "winerror"
    # attribute exclusive to such intsances.
    return exception.winerror == _ERROR_INVALID_NAME

# ....................{ TESTERS ~ version                 }....................
@func_cached
def is_version_10_or_newer() -> bool:
    '''
    ``True`` only if the current platform is Windows >= 10 (i.e., either
    Windows 10 or a newer version of Windows).
    '''

    # Avoid circular import dependencies.
    from betse.util.os import kernels

    # Return true only if the current Windows kernel version is at least 10.
    return kernels.is_version_greater_than_or_equal_to('10.0.0')

# ....................{ GETTERS                           }....................
@func_cached
def get_api_version() -> str:
    '''
    Human-readable ``.``-delimited version specifier of the Windows API
    (WinAPI, Win32) underlying the current Windows installation (e.g.,
    ``6.2.9200``).
    '''

    # If the current platform is *NOT* Windows, raise an exception.
    die_unless_windows()

    # Return the second item of the 4-tuple "(release, version, csd, ptype)"
    # returned by this function.
    return platform.win32_ver()[1]
