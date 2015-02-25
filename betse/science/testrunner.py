#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

from betse.science.runner import SimRunner

# If this module is imported from the command line, run; else, noop.
if __name__ == '__main__':
    sim_runner = SimRunner()
    sim_runner.simulate()
