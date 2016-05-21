#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Global test configuration for all tests.

`py.test` implicitly imports _all_ functionality defined by this module into
_all_ test modules.
'''

# ....................{ IMPORTS ~ fixture                  }....................
from betse_test.fixture.ignition import betse_init
