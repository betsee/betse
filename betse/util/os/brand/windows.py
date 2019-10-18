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

#FIXME: Generalize this configuration to support PowerShell. Sadly, PowerShell
#is fundamentally insane. (Who didn't see that one coming?) In particular,
#PowerShell insanely interprets *ANY* attempt to write to stderr as a fatal
#exception. See our "appveyor.yml" configuration for further details.
#
#Since Microsoft clearly has no compelling interest in resolving this blatantly
#broken behaviour, we *MUST* do so on their behalf. To do so, consider:
#
#* Implementing a new betse.util.os.brand.windows.is_shell_powershell() tester
#  returning True only if the current shell environment is PowerShell. Note
#  that, as Windows only natively supports two shell environments, the related
#  is_shell_cmd() tester is trivially implementable as follows:
#    @func_cached
#    def is_shell_cmd() -> bool:
#        return not is_shell_powershell()
#* Generalize the _init_logger_root_handler_std() method defined below to:
#
#       # Expand this unconditional assignment...
#       self._logger_root_handler_stderr = StreamHandler(sys.stderr)
#
#       # ...into these conditional assignments.
#       if windows.is_windows() and windows.is_shell_powershell():
#           self._logger_root_handler_stderr = StreamHandler(sys.stdout)
#       else:
#           self._logger_root_handler_stderr = StreamHandler(sys.stderr)
#
#In optimistic theory, the above should suffice. </apathetic_shrug>
#FIXME: Actually, to implement this sanely, we'll need to generalize our use of
#"sys.stdout" and "sys.stderr" as follows:
#
#* Define a new betse.util.io.stderrs.get_stderr() function resembling:
#    @func_cached
#    def get_stderr() -> ????:
#       return (
#           sys.stdout
#           if windows.is_windows() and windows.is_shell_powershell():
#           sys.stderr)
#* Define a new betse.util.io.stdouts.get_stdout() function resembling:
#    def get_stdout() -> ????:
#        return sys.stdout
#* Replace all references to "sys.stdout" with calls to
#  betse.util.io.stdouts.get_stdout().
#* Replace all references to "sys.stderr" with calls to
#  betse.util.io.stderrs.get_stderr().
#FIXME: Unfortunately, implementing the requisite is_shell_powershell() tester
#in this submodule will prove highly non-trivial. The reason, of course, is
#that the Windows API provides no convenient means of differentiating a child
#PowerShell from CMD.exe process. That said, doing so is absolutely feasible.
#
#The standard solution to this issue is to leverage the third-party "psutil"
#package to query whether the name of the parent process matches that of a
#known PowerShell executable, as detailed here:
#    https://stackoverflow.com/a/55598796/2809027
#
#Since that package is C-based and hence heavyweight, however, the above
#solution should *ONLY* be pursued if everything written below fails. So, what
#alternatives do we happen? Simple! We implement is_shell_powershell() in pure
#Python using only standard packages and modules. But how do we do that? Here
#we go:
#
#* Invoke the "tasklist" executable effectively guaranteed to be available
#  in the current %PATH% in all modern (i.e., post-2012) Windows installations.
#* Pass that executable options that render the output amenable to trivial
#  regular expression-based parsing ala the official "talklist" documentation:
#    https://docs.microsoft.com/en-us/windows-server/administration/windows-commands/tasklist
#* Parse the "stdout" emitted by that executable using Python's standard "csv"
#  module ala this unofficial, clever gist:
#    https://gist.github.com/intco/6149781

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
    ``True`` only if the current platform is **Cygwin Microsoft Windows**
    (i.e., running the Cygwin POSIX compatibility layer).
    '''

    return sys.platform == 'cygwin'


@func_cached
def is_windows_vanilla() -> bool:
    '''
    ``True`` only if the current platform is **vanilla Microsoft Windows**
    (i.e., *not* running the Cygwin POSIX compatibility layer).
    '''

    return sys.platform == 'win32'


@func_cached
def is_windows_wsl() -> bool:
    '''
    ``True`` only if the current platform is **Windows Subsystem for Linux
    (WSL)** (i.e., the Microsoft-flavoured Linux kernel optionally supported by
    Windows 10).

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
