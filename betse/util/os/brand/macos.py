#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2022 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Apple macOS-specific facilities.
'''

# ....................{ IMPORTS                            }....................
import platform
from betse.exceptions import BetseOSException
from betse.util.io.log.logs import log_warning
from betse.util.type.decorator.decmemo import func_cached
from ctypes import CDLL, byref, c_int

# ....................{ CONSTANTS                          }....................
_SECURITY_FRAMEWORK_DYLIB_FILENAME = (
    '/System/Library/Frameworks/Security.framework/Security')
'''
Absolute path of the system-wide ``"Security.framework"`` Macho-O shared library
providing the macOS-specific security context for the current process.

This library is dynamically loadable into the address space of the current
process with the :class:`ctypes.CDLL` class. Since all Macho-O shared libraries
necessarily have the filetype ``".dylib"``, this filetype can be safely omitted.
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

    See Also
    ----------
    :func:`is_macos`
        Further details.
    '''

    # Avoid circular import dependencies.
    from betse.util.os import oses

    # If the current platform is *NOT* macOS, raise an exception.
    if not is_macos():
        raise BetseOSException(f'{oses.get_name()} not macOS.')

# ....................{ TESTERS                            }....................
@func_cached
def is_macos() -> bool:
    '''
    ``True`` only if the current platform is Apple macOS, the operating
    system previously known as "OS X."
    '''

    return platform.system() == 'Darwin'


@func_cached
def is_aqua() -> bool:
    '''
    ``True`` only if the current process has access to the Aqua display server
    specific to macOS, implying this process to be headfull and hence support
    both CLIs and GUIs.

    See Also
    ----------
    https://developer.apple.com/library/content/technotes/tn2083/_index.html#//apple_ref/doc/uid/DTS10003794-CH1-SUBSECTION19
        "Security Context" subsection of "Technical Note TN2083: Daemons and
        Agents," a psuedo-human-readable discussion of the
        ``sessionHasGraphicAccess`` bit flag returned by the low-level
        ``SessionGetInfo()`` C function.
    '''

    # Avoid circular import dependencies.
    from betse.util.path.files import is_file
    from betse.util.os.command.cmdexit import SUCCESS

    # If the current platform is *NOT* macOS, return false.
    if not is_macos():
        return False
    # Else, the current platform is macOS.

    # If the system-wide Macho-O shared library providing the macOS
    # security context for the current process does *NOT* exist (after
    # following symbolic links)...
    if not is_file(_SECURITY_FRAMEWORK_DYLIB_FILENAME):
        # Emit a non-fatal warning. Theoretically, this shared library should
        # *ALWAYS* exist across all macOS versions (including those still
        # actively maintained as of 2022 Q2). Pragmatically, this shared library
        # appears to *NOT* exist (for unknown reasons) on GitHub Actions macOS
        # runners. He have no control over GitHub Actions. Let's complain! \o/
        log_warning(
            'macOS shared library "%s" not found.',
            _SECURITY_FRAMEWORK_DYLIB_FILENAME)

        # Return false.
        return False
    # Else, this shared library exists.

    # Attempt all of the following in a safe manner catching, logging, and
    # converting exceptions into a false return value. This tester is *NOT*
    # mission-critical and hence should *NOT* halt the application on
    # library-specific failures.
    try:
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
    # If the above logic failed with any exception...
    except Exception as exception:
        # Human-readable exception message harvested from this exception,
        # defined as either:
        # * If this exception is a platform-specific "OSError" (as is likely due
        #   to calling low-level platform-specific macOS kernel functions
        #   above), the "OSError.strerror" instance variable of this exception.
        # * Else, the standard string representation of this exception.
        exception_message = getattr(exception, 'strerror', str(exception))

        # Log a non-fatal warning informing users of this failure.
        log_warning(
            'macOS SessionGetInfo() function failed, as %s.', exception_message)

    # Assume this process to *NOT* have access to the Aqua display server.
    return False
