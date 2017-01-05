#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Global functional test configuration common to all interfaces (e.g., CLI, GUI).

`py.test` implicitly imports all functionality defined by this module into all
functional test modules. As this functionality includes all publicly declared
functional fixtures in this `fixture` subpackage, these tests may reference
these fixtures without explicit imports.
'''

# ....................{ IMPORTS ~ fixture                  }....................
