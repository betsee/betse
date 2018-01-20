#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Microsoft Windows-specific facilities.

Caveats
----------
Operating system-specific logic is poor form and should be leveraged only where
necessary.
'''

# ....................{ IMPORTS                            }....................
from betse.exceptions import BetseOSException
from betse.util.type.types import type_check, NoneType

# ....................{ TYPES                              }....................
WindowsErrorType = NoneType
'''
Windows-specific :class:`WindowsError` class if the current platform is Windows
_or_ the :class:`NoneType` class otherwise.

Since all exceptions subclass the root :class:`Exception` superclass _and_
since the :class:`NoneType` class does not do so, this class is only an
exception class under Windows. In particular, attempting to either catch this
class as an exception _or_ to test whether a caught exception is an instance of
this class is guaranteed to both be safe and behave as expected regardless of
the current platform. In short: don't worry, be emoji.
'''

# For use in type checking, this class is safely defined at the top-level in a
# platform-agnostic manner.
try:
    WindowsErrorType = WindowsError
except:
    pass

# ....................{ CONSTANTS ~ error codes            }....................
# For conformance, the names of all error code constants defined below are
# exactly as specified by Microsoft itself. Sadly, Python fails to provide these
# magic numbers for us.

_ERROR_INVALID_NAME = 123
'''
Microsoft Windows-specific error code indicating an invalid pathname.

See Also
----------
https://msdn.microsoft.com/en-us/library/windows/desktop/ms681382%28v=vs.85%29.aspx
    Official listing of all such codes.
'''

# ....................{ EXCEPTIONS                         }....................
def die_unless_windows() -> None:
    '''
    Raise an exception unless the current platform is Microsoft Windows.
    '''

    # Avoid circular import dependencies.
    from betse.util.os import oses

    # If the current platform is *NOT* Windows, raise an exception.
    if not oses.is_windows():
        raise BetseOSException(
            'Current platform {} not Windows.'.format(oses.get_name()))

# ....................{ TESTERS                            }....................
@type_check
def is_exception_pathname_invalid(exception: WindowsErrorType) -> bool:
    '''
    `True` only if the passed Windows-specific exception was raised from an
    erroneous attempt to read or write an invalid pathname.
    '''

    # Raise an exception unless the current platform is Windows.
    die_unless_windows()

    # The above type checking guarantees this exception to be an instance
    # of the Windows-specific "WindowsError" class. Test the "winerror"
    # attribute exclusive to such intsances.
    return exception.winerror == _ERROR_INVALID_NAME
