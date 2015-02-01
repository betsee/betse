#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''Main entry point of `betse`'s command line interface (CLI).

This module is a thin wrapper intended to be imported only from the command line
(e.g., via `python -m betse`). This module is *not* imported by `setuptools`'s
`entry_points` facility and hence is *not* imported by external scripts
managed by `setuptools` (e.g., `betse`).
'''

# ....................{ IMPORTS                            }....................
from betse.cli.cli import main

# ....................{ MAIN                               }....................
# If this module is imported from the command line, run betse's CLI; else, noop.
if __name__ == '__main__':
    main()
