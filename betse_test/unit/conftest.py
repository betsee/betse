#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Global unit test configuration for BETSE.

:mod:`pytest` implicitly imports all functionality defined by this module into
all unit test modules. As this functionality includes all publicly declared
functional fixtures in this fixture subpackage, these tests may reference these
fixtures without explicit imports.
'''

# ....................{ IMPORTS ~ fixture : autouse       }....................
# Import fixtures automatically run at the start of the current test session,
# typically *NOT* manually required by specific tests, *AFTER* importing all
# non-autouse fixtures possibly required by these autouse fixtures above.

from betse_test.fixture.initter import betse_init_package
