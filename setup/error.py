#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''Error-handling functions for `betse`-specific `setuptools` commands.'''

# ....................{ IMPORTS                            }....................
from distutils.errors import DistutilsPlatformError
import os

# ....................{ EXCEPTIONS                         }....................
def raise_exception_if_os_non_posix():
    '''
    Raise a fatal exception if the current operating system does `not` comply
    with POSIX standards (e.g., as required for symbolic link manipulation).

    Typically, this implies such system to be Windows.
    '''
    if os.name != 'posix':
        raise DistutilsPlatformError(
            'This command requires POSIX compliance. Distressingly, the current '
            'operating system is POSIX-noncompliant (e.g., Windows).'
        )

# --------------------( WASTELANDS                         )--------------------
