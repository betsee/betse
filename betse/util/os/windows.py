#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Microsoft Windows-specific facilities.

Caveats
----------
Operating system-specific logic is poor form and should be leveraged only where
necessary.
'''

# ....................{ IMPORTS                            }....................

# ....................{ CONSTANTS ~ error codes            }....................
# For conformance, the names of all error code constants defined below are
# exactly as specified by Microsoft itself. Sadly, Python fails to provide these
# magic numbers for us.

ERROR_INVALID_NAME = 123
'''
Microsoft Windows-specific error code indicating an invalid pathname.

See Also
----------
https://msdn.microsoft.com/en-us/library/windows/desktop/ms681382%28v=vs.85%29.aspx
    Official listing of all such codes.
'''
