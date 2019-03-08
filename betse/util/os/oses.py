#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level operating system (OS) facilities.

Caveats
----------
**Operating system-specific logic is poor form.** Do so *only* where necessary.
'''

# ....................{ IMPORTS                           }....................
import os, platform, sys
from betse import metadata
from betse.util.io.log import logs
from betse.util.type.decorator.decmemo import func_cached
from betse.util.type.iterable.mapping.mapcls import OrderedArgsDict

# ....................{ INITIALIZERS                      }....................
def init() -> None:
    '''
    Validate the current operating system.

    This function (in order):

    #. Logs a non-fatal warning if this operating system is a non-WSL variant
       of Microsoft Windows (e.g., vanilla Windows, Cygwin Windows).
    #. Logs a non-fatal warning if this operating system is *not* recognized as
       officially supported by this application (e.g., BSD*, Solaris).
    '''

    # Log this validation.
    logs.log_debug('Validating operating system...')

    # Human-readable string describing the set of all officially supported
    # platforms known to interoperate sanely with this application.
    supported_oses = (
        'Consider running {name} only under Linux or macOS. '
        'Note that Linux is now freely emulatable under Windows 10 '
        'via the Windows Subsystem for Linux (WSL). '
        'See official installation instructions at:\n'
        '\thttps://msdn.microsoft.com/en-us/commandline/wsl/install-win10'.format(
            name=metadata.NAME))

    # If this is a non-WSL Windows variant, log a non-fatal warning.
    if is_windows():
        logs.log_warning(
            'Windows operating system detected. '
            'Python itself and third-party scientific frameworks for Python '
            '(e.g., Numpy, Matplotlib) are well-known to behave suboptimally '
            'under Windows, '
            'impeding the reliability and scalability of modelling. %s',
            supported_oses)

    # If this platform is officially unsupported by this application, log a
    # non-fatal warning.
    if not is_supported():
        logs.log_warning(
            'Unsupported operating system "%s" detected. %s',
            get_name(), supported_oses)

# ....................{ TESTERS                           }....................
@func_cached
def is_supported() -> bool:
    '''
    ``True`` only if the current operating system is officially supported by
    this application.

    This function currently only returns ``True`` for the following platforms
    (in no particular order):

    * Apple macOS.
    * Linux.
    * Microsoft Windows.

    Caveats
    ----------
    This function returning ``True`` does *not* necessarily imply this
    application to behave optimally under this operating system. In particular,
    Microsoft Windows is officially supported by popular demand but well-known
    to behave suboptimally with respect to Python itself and third-party
    scientific frameworks for Python (e.g., Numpy, Matplotlib).
    '''

    return is_linux() or is_macos() or is_windows()

# ....................{ TESTERS ~ posix                   }....................
@func_cached
def is_posix() -> bool:
    '''
    ``True`` only if the current operating system complies with POSIX standards
    (e.g., as required for POSIX-compliant symbolic link support).

    Typically, this implies this system to *not* be vanilla Microsoft Windows
    and hence to be either:

    * A genuinely POSIX-compliant system.
    * A Cygwin-based Windows application (e.g., CLI terminal, GUI application).
    '''

    return os.name == 'posix'


@func_cached
def is_linux() -> bool:
    '''
    ``True`` only if the current operating system is either Linux or an
    operating system successfully masquerading to a sufficiently accurate
    degree as Linux (e.g., the Windows Subsystem for Linux (WSL) but *not*
    Cygwin Windows).
    '''

    return platform.system() == 'Linux'


@func_cached
def is_macos() -> bool:
    '''
    ``True`` only if the current operating system is Apple macOS, the operating
    system previously known as "OS X."
    '''

    return platform.system() == 'Darwin'

# ....................{ TESTERS ~ windows                 }....................
@func_cached
def is_windows() -> bool:
    '''
    ``True`` only if the current operating system is Microsoft Windows.

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
    ``True`` only if the current operating system is **Cygwin Microsoft
    Windows** (i.e., running the Cygwin POSIX compatibility layer).
    '''

    return sys.platform == 'cygwin'


@func_cached
def is_windows_vanilla() -> bool:
    '''
    ``True`` only if the current operating system is **vanilla Microsoft
    Windows** (i.e., *not* running the Cygwin POSIX compatibility layer).
    '''

    return sys.platform == 'win32'


@func_cached
def is_windows_wsl() -> bool:
    '''
    ``True`` only if the current operating system is **Windows Subsystem for
    Linux (WSL)** (i.e., the Microsoft-flavoured Linux kernel optionally
    supported by Windows 10).

    See Also
    ----------
    https://www.reddit.com/r/bashonubuntuonwindows/comments/85jghk/how_to_allow_python_to_know_im_on_windows/dvxwy9t
        Reddit post strongly inspiring this implementation.
    '''

    # If this the active Python interpreter is *NOT* operating under a Linux
    # kernel, return false immediately.
    if not is_linux():
        return False

    # Flavour of this Linux kernel.
    kernel_flavour = platform.uname()[3]

    # Return true only if this is a Microsoft-flavoured Linux kernel.
    return 'microsoft' in kernel_flavour

# ....................{ GETTERS                           }....................
@func_cached
def get_name() -> str:
    '''
    Human-readable name of the current operating system.

    This function returns:

    * Under Linux, the name of this Linux distribution, detected as follows:

      * If the :func:`platform.linux_distribution` function is available, the
        first element of the 3-tuple returned by this function (e.g.,
        ``CentOS``).
      * Else, the string returned by the :func:`platform.system` function.
        Since this is probably the generic string ``Linux``, this is only a
        fallback.

    * Under macOS, ``macOS``.
    * Under Windows, the flavour of the operating system-level kernel exposed
      by the shell environment to which the active Python interpreter is
      attached. Specifically, if this is:

      * The Windows Subsystem for Linux (WSL), ``Windows (WSL)``. Note,
        however, that the WSL masquerades as Linux. Thankfully, it does so to a
        remarkably more reliable degree than Cygwin Windows.
      * Cygwin Windows, ``Windows (Cygwin)``. Note, however, that Cygwin
        masquerades as Linux, replicating the POSIX process model and GNU/Linux
        userland to a remarkable degree -- subject to the Big List of Dodgy
        Apps (BLODA), proprietary limitations of the Win32 API, and similar
        caveats.
      * Under vanilla Windows shell environments (e.g., CMD.exe, PowerShell),
        merely ``Windows``.

    * Under all other platforms, the string returned by the
      :func:`platform.system` function.
    '''

    # Name to be returned, defaulting to that returned by platform.system().
    # Since this typically corresponds to the low-level name of the current
    # kernel (e.g., "Darwin", "Linux") rather than the high-level name of the
    # current OS (e.g., "macOS", "CentOS"), this is only a fallback.
    os_name = platform.system()

    # If Linux...
    if is_linux():
        # If the Linux-specific platform.linux_distribution() function is
        # available, return the first item of the 3-tuple
        # (distname, version, id) returned by this function (e.g., "CentOS").
        # Since platform.system() returns the low-level kernel name "Linux",
        # this name is ignored where feasible.
        if hasattr(platform, 'linux_distribution'):
            os_name = platform.linux_distribution()[0]
        # Else if this is actually the Windows Subsystem for Linux (WSL)
        # masquerading as Linux, return "Windows (WSL)" rather than "Linux".
        elif is_windows_wsl():
            os_name = 'Windows (WSL)'
        # Else, reuse the name returned by the prior call to platform.system().
    # If macOS, return "macOS". Since platform.system() returns the low-level
    # kernel name "Darwin", this name is ignored.
    elif is_macos():
        os_name = 'macOS'
    # If Cygwin Windows, return "Windows (Cygwin)". Since platform.system()
    # returns a non-human-readable low-level uppercase label specific to the
    # word size of the current Python interpreter (e.g., "CYGWIN_NT-5.1*"),
    # this label is ignored.
    elif is_windows_cygwin():
        os_name = 'Windows (Cygwin)'
    # If non-Cygwin Windows, return merely "Windows". The name returned by
    # platform.system() depends on the current Windows version as follows:
    #
    # * If this is Windows Vista, this name is "Microsoft".
    # * Else, this name is "Windows".
    #
    # Hence, this name is ignored.
    elif is_windows_vanilla():
        os_name = 'Windows'
    # Else, reuse the name returned by the prior call to platform.system().

    # Return the name established above.
    return os_name


@func_cached
def get_version() -> str:
    '''
    Human-readable ``.``-delimited version specifier string of the current
    operating system.

    This function returns:

    * Under Linux, the version reported by this Linux distribution as follows:

      * If the :func:`platform.linux_distribution` function is available, the
        second element of the 3-tuple returned by this function (e.g.,
        ``6.4``).
      * Else, the string returned by the :func:`platform.release` function.
        Since this is probably the non-human-readable ``.``- and
        ``-``-delimited version specifier string of the current Linux kernel
        (e.g., ``4.1.15-gentoo-r1``), this is only a fallback.

    * Under macOS, this system's current major and minor version (e.g.,
      ``10.9``).
    * Under Windows, this system's current major version (e.g., ``8``).
      Unfortunately, there appears to be no consistently reliable means of
      obtaining this system's current major *and* minor version (e.g.,
      ``8.1``).
    * Under all other platforms, the string returned by the
      :func:`platform.release` function.
    '''

    # Version specifier to be returned, defaulting to that returned by
    # platform.release(). Since this typically corresponds to the low-level
    # version of the current kernel (e.g., "4.1.15") rather than the high-level
    # version of the current OS (e.g., "6.4"), this is only a fallback.
    os_version = platform.release()

    # If the Linux-specific platform.linux_distribution() function is
    # available, return the second element of the 3-tuple "(distname, version,
    # id)" returned by this function (e.g., "6.4"). Since platform.release()
    # returns the low-level kernel version (e.g., "4.1.15"), this version is
    # ignored when feasible.
    if is_linux() and hasattr(platform, 'linux_distribution'):
        os_version = platform.linux_distribution()[1]
    # If macOS *AND* the platform.mac_ver() function is available, return the
    # first element of the 3-tuple returned by this function (e.g., "10.9").
    # Since platform.release() returns the low-level kernel version
    # (e.g., "13.0.0"), this version is ignored when feasible.
    elif is_macos() and hasattr(platform, 'mac_ver'):
        # platform.mac_ver() returns a 3-tuple (release, versioninfo, machine)
        # of macOS-specific metadata, where "versioninfo" is itself a 3-tuple
        # (version, dev_stage, non_release_version). Return the high-level
        # "release" element (e.g., "10.9") rather than the optionally defined
        # low-level "version" element of the "versioninfo" element, which is
        # typically the empty string and hence useless.
        os_version = platform.mac_ver()[0]
    # If Windows, accept the default version specifier returned by
    # platform.release() (e.g., "8", "10"). Since this specifier appears to
    # coincide exactly with the first element of the 4-tuple returned by the
    # conditionally available platform.win32_ver() function, there is no
    # benefit to the latter approach.

    # Return this version.
    return os_version

# ....................{ GETTERS ~ metadata                }....................
def get_metadata() -> OrderedArgsDict:
    '''
    Ordered dictionary synopsizing the current operating system.
    '''

    # Return this dictionary.
    return OrderedArgsDict(
        'name', get_name(),
        'version', get_version(),
    )
