#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Apple macOS-specific facilities.
'''

# ....................{ IMPORTS                            }....................
from betse.exceptions import BetseOSException
from betse.util.io.log import logs
from ctypes import CDLL, byref, c_int

# ....................{ CONSTANTS                          }....................
_SECURITY_FRAMEWORK_DYLIB_FILENAME = (
    '/System/Library/Frameworks/Security.framework/Security')
'''
Absolute path of the system-wide `Security.framework` Macho-O shared library
providing the macOS-specific security context for the current process.

This library is dynamically loadable into the address space of the current
process with the :class:`ctypes.CDLL` class. Since all Macho-O shared libraries
necessarily have the filetype `dylib`, this filetype is safely omitted here.
'''


_SECURITY_SESSION_ID_CURRENT = -1
'''
Magic integer defined as `callerSecuritySession` by the macOS-specific
`/System/Library/Frameworks/Security.Framework/Headers/AuthSession.h` C header
suitable for passing to C functions accepting parameters of C type
`SecuritySessionId` (e.g., `SessionGetInfo()`).

When passed, this integer signifies the **current security session** (i.e., the
the security session to which the current process belongs).

See Also
----------
https://opensource.apple.com/source/libsecurity_authorization/libsecurity_authorization-32564/lib/AuthSession.h
    C header defining this magic integer.
'''


_SECURITY_SESSION_HAS_GRAPHIC_ACCESS = 0x0010
'''
Bit flag defined as `sessionHasGraphicAccess` by the macOS-specific
`/System/Library/Frameworks/Security.Framework/Headers/AuthSession.h` C header
masking the attributes bit field returned by the `SessionGetInfo()` C function
also declared by that header.

When enabled, this bit signifies the current process to have access to the Aqua
display server and hence be headfull (rather than headless).

See Also
----------
https://opensource.apple.com/source/libsecurity_authorization/libsecurity_authorization-32564/lib/AuthSession.h
    C header defining this bit flag.
'''

# ....................{ EXCEPTIONS                         }....................
def die_unless_macos() -> None:
    '''
    Raise an exception unless the current platform is Apple macOS.
    '''

    # Avoid circular import dependencies.
    from betse.util.os import oses

    # If the current platform is *NOT* macOS, raise an exception.
    if not oses.is_macos():
        raise BetseOSException(
            'Current platform {} not macOS.'.format(oses.get_name()))

# ....................{ TESTERS                            }....................
def is_aqua() -> bool:
    '''
    `True` only if the current process has access to the Aqua display server
    specific to macOS, implying this process to be headfull and hence support
    both CLIs and GUIs.

    See Also
    ----------
    https://developer.apple.com/library/content/technotes/tn2083/_index.html#//apple_ref/doc/uid/DTS10003794-CH1-SUBSECTION19
        "Security Context" subsection of "Technical Note TN2083: Daemons and
        Agents," a psuedo-human-readable discussion of the
        `sessionHasGraphicAccess` bit flag returned by the low-level
        `SessionGetInfo()` C function.
    '''

    # Avoid circular import dependencies.
    from betse.util.path import files
    from betse.util.path.command.cmdexit import SUCCESS

    # Raise an exception unless the current platform is macOS.
    die_unless_macos()

    # Attempt all of the following in a safe manner catching, logging, and
    # converting exceptions into a false return value. This tester is *NOT*
    # mission-critical and hence should *NOT* halt the application on
    # library-specific failures.
    try:
        # If the system-wide Macho-O shared library providing the macOS
        # security context for the current process does *NOT* exist (after
        # following symbolic links), raise an exception.
        files.die_unless_file(_SECURITY_FRAMEWORK_DYLIB_FILENAME)

        # Dynamically load this library into the address space of this process.
        security_framework = CDLL(_SECURITY_FRAMEWORK_DYLIB_FILENAME)

        # Possibly non-unique identifier of the security session to request the
        # attributes of, signifying that of the current process.
        session_id = _SECURITY_SESSION_ID_CURRENT

        # Unique identifier of the requested security session, returned
        # by reference from the SessionGetInfo() C function called below. This
        # identifier is useless for our purposes and hence ignored below.
        session_id_real = c_int(0)

        # Attributes bit field of the requested security session, returned by
        # reference from the SessionGetInfo() C function called below.
        session_attributes = c_int(0)

        # C-style error integer returned by calling the SessionGetInfo() C
        # function exported by this Macho-O shared library, passing:
        #
        # * The input non-unique session identifier by value.
        # * The output unique session identifier by reference.
        # * The output session attributes integer by reference.
        session_errno = security_framework.SessionGetInfo(
            session_id, byref(session_id_real), byref(session_attributes))

        # This process has access to the Aqua display server if and only if...
        return (
            # The above function call succeeded *AND*...
            session_errno == SUCCESS and
            # The session attributes bit field returned by this call has the
            # corresponding bit flag enabled.
            session_attributes.value & _SECURITY_SESSION_HAS_GRAPHIC_ACCESS
        )

    # If the above logic fails with any exception...
    except Exception as exc:
        # Log a non-fatal warning informing users of this failure.
        logs.log_warning(
            'macOS-specific SessionGetInfo() C function failed: {}'.format(
                exc.strerror))

        # Assume this process to *NOT* have access to the Aqua display server.
        return False
