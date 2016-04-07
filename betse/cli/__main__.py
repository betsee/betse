#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Main entry point of `betse`'s command line interface (CLI).

This module is a thin wrapper intended to be imported only from the command line
(e.g., via `python -m betse`). This module is *not* imported by `setuptools`'s
`entry_points` facility and hence is *not* imported by external scripts
managed by `setuptools` (e.g., `betse`).
'''

# ....................{ IMPORTS                            }....................
import sys
from betse.cli.clicli import CLICLI

# ....................{ MAIN                               }....................
def main() -> int:
    '''Run `betse`'s command line interface (CLI).

    This function is provided as a convenience to callers requiring procedural
    functions rather than conventional methods (e.g., `setuptools`).

    Returns
    ----------
    int
        Exit status of such interface. This is a non-negative integer in
        `[0, 255]` where 0 signifies success and all other values failure.
    '''
    # print('In main')
    return CLICLI().run()

# ....................{ MAIN                               }....................
# If this module is imported from the command line, run BETSE's CLI; else, noop.
#
# For POSIX compliance, the value returned by this function (ideally a single-
# byte integer) will be propagated back to the calling shell as this script's
# exit status.
if __name__ == '__main__':
    sys.exit(main())
