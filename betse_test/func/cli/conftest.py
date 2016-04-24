#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Global functional test configuration for BETSE's command-line interface (CLI).

`py.test` implicitly imports all functionality defined by this module into all
functional test modules. As this functionality includes all publicly declared
functional fixtures in the `fixture` subpackage, functional tests may depend on
functional fixtures without explicitly importing those fixtures.
'''

# ....................{ IMPORTS ~ fixture                  }....................
# from betse_test_func.util.context.context import SimTestContext
