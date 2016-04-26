#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Main entry point of `betse`'s command line interface (CLI).

This module is a thin wrapper intended to be imported only from the command line
(e.g., via `python -m betse`). This module is *not* imported by `setuptools`'s
`entry_points` facility and hence is *not* imported by external scripts
managed by `setuptools` (e.g., `betse`).
'''

# ....................{ IMPORTS                            }....................
from betse.cli.clicli import CLICLI
from betse.util.command import exits

# ....................{ MAIN                               }....................
def main(arg_list: list = None) -> int:
    '''
    Run BETSE's command-line interface (CLI) with the passed arguments if
    non-`None` _or_ with the arguments passed on the command line (i.e.,
    `sys.argv`) otherwise.

    This function is provided as a convenience to callers requiring procedural
    functions rather than conventional methods (e.g., `setuptools`).

    Parameters
    ----------
    arg_list : list
        List of zero or more arguments to pass to this interface. Defaults to
        `None`, in which case arguments passed on the command line (i.e.,
        `sys.argv`) will be used instead.

    Returns
    ----------
    int
        Exit status of this interface, in the range `[0, 255]`.
    '''

    return CLICLI().run(arg_list)

# ....................{ MAIN                               }....................
# If this module is imported from the command line, run BETSE's CLI; else, noop.
# For POSIX compliance, the exit status returned by this function is propagated
# to the caller as this script's exit status.
if __name__ == '__main__':
    exits.exit(main())
