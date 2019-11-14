#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Subsidiary entry point of this application's command line interface (CLI),
preserving backward compatibility with prior versions.
'''

# ....................{ IMPORTS                           }....................
from betse.__main__ import main
from betse.util.os.command import cmdexit

# ....................{ MAIN                              }....................
# If this module is imported from the command line, run this application's CLI;
# else, noop. For POSIX compliance, the exit status returned by this function
# is propagated to the caller as this script's exit status.
if __name__ == '__main__':
    cmdexit.exit_with_status(main())
