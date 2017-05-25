#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level shell-specific string facilities.
'''

# ....................{ IMPORTS                            }....................
from betse.util.type.types import type_check

# ....................{ QUOTERS                            }....................
@type_check
def shell_quote(text: str) -> str:
    '''
    Shell-quote the passed string in a platform-specific manner.

    If the current platform is:

    * *Not* Windows (e.g., Linux, OS X), the returned string is guaranteed to be
      suitable for passing as an arbitrary positional argument to external
      commands.
    * Windows, the returned string is suitable for passing *only* to external
      commands parsing arguments in the same manner as the Microsoft C runtime.
      While *all* applications on POSIX-compliant systems are required to parse
      arguments in the same manner (i.e., according to Bourne shell lexing), no
      such standard applies to Windows applications. Shell quoting is therefore
      fragile under Windows -- like pretty much everything.
    '''

    # Avoid circular import dependencies.
    from betse.util.os import oses

    # If the current OS is Windows, do *NOT* perform POSIX-compatible quoting.
    # Windows is POSIX-incompatible and hence does *NOT* parse command-line
    # arguments according to POSIX standards. In particular, Windows does *NOT*
    # treat single-quoted arguments as single arguments but rather as multiple
    # shell words delimited by the raw literal `'`. This is circumventable by
    # calling an officially undocumented Windows-specific function. (Awesome.)
    if oses.is_windows():
        import subprocess
        return subprocess.list2cmdline([text])
    # Else, perform POSIX-compatible quoting.
    else:
        import shlex
        return shlex.quote(text)
