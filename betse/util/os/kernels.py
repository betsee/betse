#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level **kernel** (i.e., central core software of the current operating
system) facilities.

Caveats
----------
**Operating system-specific logic is poor form.** Do so *only* where necessary.
'''

# ....................{ IMPORTS                           }....................
import platform
from betse.util.type.decorator.decmemo import func_cached
from betse.util.type.iterable.mapping.mapcls import OrderedArgsDict
from betse.util.type.types import type_check, VersionTypes

# ....................{ TESTERS                           }....................
@type_check
def is_version_greater_than_or_equal_to(version: VersionTypes) -> bool:
    '''
    ``True`` only if the current kernel version is greater than or equal to the
    passed version.

    Parameters
    ----------
    version : VersionTypes
        Version (e.g., ``1.0.2``, ``(1, 0, 2)``) to test the current kernel
        version against.

    Returns
    ----------
    bool
        ``True`` only if the current kernel version is greater than or equal to
        this version.

    Raises
    ----------
    pkg_resources.packaging.version.InvalidVersion
        If this version is *not* `PEP 440-compliant`_.

    .. _PEP 440-compliant:
       https://www.python.org/dev/peps/pep-0440

    See Also
    ----------
    :func:`get_version`
        Further details.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.numeric import versions

    # Current kernel version as a human-readable string.
    kernel_version = get_version()

    # Return true only if this version is greater than or equal to the passed
    # version.
    return versions.is_greater_than_or_equal_to(kernel_version, version)

# ....................{ GETTERS                           }....................
@func_cached
def get_name() -> str:
    '''
    Machine-readable name of the current kernel.

    While preserving human-readability is helpful, machine-readability is the
    principal use case for this function's return value (e.g., for use in
    dictionaries keyed on kernel type). This function returns:

    * Under Linux, ``Linux``.
    * Under macOS, ``Darwin``.
    * Under Windows (both Cygwin and vanilla), ``Windows``.
    * Under all other platforms, the string returned by the
      :func:`platform.system` function.
    '''

    # Return platform.system() as is, which appears to exactly corresponding to
    # the name of the current kernel on all platforms (e.g., "Darwin").
    return platform.system()


@func_cached
def get_version() -> str:
    '''
    Human-readable ``.``-delimited version specifier of the current kernel.

    This function returns:

    * Under Linux, the current Linux kernel version (e.g.,
      ``4.1.15-gentoo-r1``),
    * Under macOS, the current XNU kernel version (e.g., ``13.0.0``).
    * Under Windows, the current Windows kernel version (e.g., ``10.0.10240``).
    * Under all other platforms, the string returned by the
      :func:`platform.release` function.
    '''

    # Avoid circular import dependencies.
    from betse.util.os.brand import windows

    # Kernel version specifier to be returned.
    kernel_version = None

    # If the current platform is Windows, defer to the platform.version()
    # function. For only this platform, this function returns the desired
    # fine-grained Windows kernel version:
    #
    #     # Under Windows 10...
    #     >>> import platform
    #     >>> platform.release()
    #     "10"
    #     >>> platform.version()
    #     "10.0.10240"
    #     >>> platform.win32_ver()[1]
    #     "6.2.9200"
    #
    # However, note that the above does *NOT* generalize to any other platform:
    #
    #     # Under Gentoo Linux...
    #     >>> import platform
    #     >>> platform.release()
    #     "4.19.27-gentoo-r1"
    #     >>> platform.version()
    #     "#1 SMP Fri Apr 26 20:12:45 EDT 2019"
    #
    # Version specifier to be returned, defaulting to that returned by
    # platform.release(). While this typically corresponds to the low-level
    # version of the current kernel (e.g., "4.1.15" under Linux, "13.0.0" under
    # macOS), this is *NOT* necessarily the case (e.g., "8" rather than
    # "6.2.9200" under Windows). Hence, this is only a fallback.
    if windows.is_windows():
        kernel_version = platform.version()
    # Else, the current platform is *NOT* Windows. In this case, prefer the
    # version specifier returned by the platform.release() function to the
    # irrelevant timestamp returned by the platform.version() function.
    # Ironically, the platform.version() function does *NOT* actually return a
    # version specifier under most platforms. (Shaking my head.)
    else:
        kernel_version = platform.release()

    # Return this version.
    return kernel_version

# ....................{ GETTERS ~ metadata                }....................
def get_metadata() -> OrderedArgsDict:
    '''
    Ordered dictionary synopsizing the current kernel.
    '''

    # Return this dictionary.
    return OrderedArgsDict(
        'name', get_name(),
        'version', get_version(),
    )
