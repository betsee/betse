#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2025 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Main entry point of this application's command line interface (CLI).

This submodule is a thin wrapper intended to be:

* Indirectly imported and run from external entry point scripts installed by
  setuptools (e.g., the ``betse`` command).
* Directly imported and run from the command line (e.g., via
  ``python -m betse.cli``).
'''

# ....................{ IMPORTS                            }....................
from beartype.typing import Collection
from betse.util.os.command import cmdexit

# ....................{ MAIN                               }....................
def main(arg_list: Collection[str] = None) -> int:
    '''
    Run this application's command-line interface (CLI) with the passed
    arguments if non-``None`` *or* with the arguments passed on the command
    line (i.e., :attr:`sys.argv`) otherwise.

    This function is provided as a convenience to callers requiring procedural
    functions rather than conventional methods (e.g., :mod:`setuptools`).

    Parameters
    ----------
    arg_list : Collection[str]
        Collection of zero or more string arguments to pass to this interface.
        Defaults to :data:`None`, in which case arguments passed on the command
        line (i.e., :attr:`sys.argv`) will be used instead.

    Returns
    -------
    int
        Exit status of this interface and hence this process as an unsigned
        byte (i.e., integer in the range ``[0, 255]``).
    '''

    # Defer function-specific imports.
    from betse.appmeta import BetseAppMeta
    from betse.cli.climain import BetseCLI

    # Initialize this application (excluding mandatory dependencies of this
    # application, which will be initialized later in the startup process) by
    # instantiating the BETSE-specific application metadata singleton.
    BetseAppMeta()

    # Implement this application by instantiating and running the
    # BETSE-specific command-line interface (CLI).
    return BetseCLI().run(arg_list)

# ....................{ MAIN                              }....................
# If this module is imported from the command line, run this application's CLI;
# else, noop. For POSIX compliance, the exit status returned by this function
# is propagated to the caller as this script's exit status.
if __name__ == '__main__':
    cmdexit.exit_with_status(main())
