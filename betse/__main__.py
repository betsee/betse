#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Main entry point of BETSE's command line interface (CLI).

This submodule is a thin wrapper intended to be:

* Indirectly imported and run from external entry point scripts installed by
  setuptools (e.g., the ``betse`` command).
* Directly imported and run from the command line (e.g., via
  ``python -m betse.cli``).
'''

# ....................{ IMPORTS                            }....................
from betse.cli.climain import CLIMain
from betse.util.path.command import exits

# ....................{ MAIN                               }....................
def main(arg_list: list = None) -> int:
    '''
    Run BETSE's command-line interface (CLI) with the passed arguments if
    non-``None`` *or* with the arguments passed on the command line (i.e.,
    :attr:`sys.argv`) otherwise.

    This function is provided as a convenience to callers requiring procedural
    functions rather than conventional methods (e.g., :mod:`setuptools`).

    Parameters
    ----------
    arg_list : list
        List of zero or more arguments to pass to this interface. Defaults to
        ``None``, in which case arguments passed on the command line (i.e.,
        :attr:`sys.argv`) will be used instead.

    Returns
    ----------
    int
        Exit status of this interface, in the range ``[0, 255]``.
    '''

    # print('BETSE arg list (in main): {}'.format(arg_list))
    return CLIMain().run(arg_list)

# ....................{ MAIN                               }....................
# If this module is imported from the command line, run BETSE's CLI; else, noop.
# For POSIX compliance, the exit status returned by this function is propagated
# to the caller as this script's exit status.
if __name__ == '__main__':
    exits.exit(main())
