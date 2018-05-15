#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level operating system (OS) facilities.

Caveats
----------
**Operating system-specific logic is poor form.** Do so *only* where necessary.
'''

# ....................{ IMPORTS                            }....................
import os, platform, sys
from betse.util.type.decorator.decmemo import func_cached
from betse.util.type.mapping.mapcls import OrderedArgsDict

# ....................{ TESTERS ~ posix                    }....................
@func_cached
def is_posix() -> bool:
    '''
    ``True`` only if the current operating system complies with POSIX standards
    (e.g., as required for POSIX-compliant symbolic link support).

    Typically, this implies this system to _not_ be vanilla Microsoft Windows
    and hence to be either:

    * A genuinely POSIX-compliant system.
    * A Cygwin-based Windows application (e.g., CLI terminal, GUI application).
    '''

    return os.name == 'posix'


@func_cached
def is_linux() -> bool:
    '''
    ``True`` only if the current operating system is Linux.
    '''

    return platform.system() == 'Linux'


@func_cached
def is_macos() -> bool:
    '''
    ``True`` only if the current operating system is Apple macOS, the operating
    system previously known as "macOS."
    '''

    return platform.system() == 'Darwin'

# ....................{ TESTERS ~ windows                  }....................
@func_cached
def is_windows() -> bool:
    '''
    ``True`` only if the current operating system is Microsoft Windows.

    This function reports `True` for both vanilla and Cygwin Microsoft Windows.
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

# ....................{ GETTERS                            }....................
@func_cached
def get_name() -> str:
    '''
    Human-readable name of the current operating system.

    This function returns:

    * Under Linux, the name of this Linux distribution, detected as follows:
      * If the `platform.linux_distribution()` function is available, the first
        element of the 3-tuple returned by this function (e.g., `CentOS`).
      * Else, the string returned by the `platform.system()` function. Since
        this is probably the generic string `Linux`, this is only a fallback.
    * Under macOS, `macOS`.
    * Under Cygwin Windows, `Windows (Cygwin)`. Note, however, that Cygwin
      masquerades as Linux, replicating the POSIX process model and GNU/Linux
      userland to a remarkable degree -- subject to the Big List of Dodgy Apps
      (BLODA), proprietary limitations of the Win32 API, and similar caveats.
    * Under non-Cygwin Windows, merely `Windows`.
    * Under all other platforms, the string returned by the `platform.system()`
      function.
    '''

    # Name to be returned, defaulting to that returned by platform.system().
    # Since this typically corresponds to the low-level name of the current
    # kernel (e.g., "Darwin", "Linux") rather than the high-level name of the
    # current OS (e.g., "macOS", "CentOS"), this is only a fallback.
    os_name = platform.system()

    # If the Linux-specific platform.linux_distribution() function is available,
    # return the first element of the 3-tuple "(distname, version, id)" returned
    # by this function (e.g., "CentOS"). Since platform.system() returns the
    # low-level kernel name "Linux", this name is ignored when feasible.
    if is_linux() and hasattr(platform, 'linux_distribution'):
        os_name = platform.linux_distribution()[0]
    # If macOS, return "macOS". Since platform.system() returns the low-level
    # kernel name "Darwin", this name is ignored.
    elif is_macos():
        os_name = 'macOS'
    # If Cygwin Windows, return "Windows (Cygwin)". Since platform.system()
    # returns a non-human-readable low-level uppercase label specific to the
    # word size of the current Python interpreter (e.g., "CYGWIN_NT-5.1*"), this
    # label is ignored.
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

    # Return this name.
    return os_name


@func_cached
def get_version() -> str:
    '''
    Human-readable `.`-delimited version specifier string of the current
    operating system.

    This function returns:

    * Under Linux, the version reported by this Linux distribution as follows:
      * If the `platform.linux_distribution()` function is available, the second
        element of the 3-tuple returned by this function (e.g., `6.4`).
      * Else, the string returned by the `platform.release()` function. Since
        this is probably the non-human-readable `.`- and `-`-delimited version
        specifier string of the current Linux kernel (e.g., `4.1.15-gentoo-r1`),
        this is only a fallback.
    * Under macOS, this system's current major and minor version (e.g., `10.9`).
    * Under Windows, this system's current major version (e.g., `8`).
      Unfortunately, there appears to be no consistently reliable means of
      obtaining this system's current major _and_ minor version (e.g., `8.1`).
    * Under all other platforms, the string returned by the `platform.release()`
      function.
    '''

    # Version specifier to be returned, defaulting to that returned by
    # platform.release(). Since this typically corresponds to the low-level
    # version of the current kernel (e.g., "4.1.15") rather than the high-level
    # version of the current OS (e.g., "6.4"), this is only a fallback.
    os_version = platform.release()

    # If the Linux-specific platform.linux_distribution() function is available,
    # return the second element of the 3-tuple "(distname, version, id)"
    # returned by this function (e.g., "6.4"). Since platform.release() returns
    # the low-level kernel version (e.g., "4.1.15"), this version is ignored
    # when feasible.
    if is_linux() and hasattr(platform, 'linux_distribution'):
        os_version = platform.linux_distribution()[1]
    # If macOS *AND* the platform.mac_ver() function is available, return the
    # first element of the 3-tuple returned by this function (e.g., "10.9").
    # Since platform.release() returns the low-level kernel version
    # (e.g., "13.0.0"), this version is ignored when feasible.
    elif is_macos() and hasattr(platform, 'mac_ver'):
        # platform.mac_ver() returns a 3-tuple "(release, versioninfo, machine)"
        # of macOS-specific metadata, where "versioninfo" is itself a 3-tuple
        # "(version, dev_stage, non_release_version)". Return the high-level
        # "release" element (e.g., "10.9") rather than the optionally defined
        # low-level "version" element of the "versioninfo" element, which is
        # typically the empty string and hence useless.
        os_version = platform.mac_ver()[0]
    # If Windows, accept the default version specifier returned by
    # platform.release() (e.g., "8", "10"). Since this specifier appears to
    # coincide exactly with the first element of the 4-tuple returned by the
    # conditionally available platform.win32_ver() function, there is no benefit
    # to the latter approach.

    # Return this version.
    return os_version

# ....................{ GETTERS ~ metadata                 }....................
def get_metadata() -> OrderedArgsDict:
    '''
    Ordered dictionary synopsizing the current operating system.
    '''

    # Return this dictionary.
    return OrderedArgsDict(
        'name', get_name(),
        'version', get_version(),
    )
