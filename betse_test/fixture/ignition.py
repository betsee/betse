#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
BETSE initialization fixtures.
'''

# ....................{ IMPORTS                            }....................
from pytest import fixture

# ....................{ FIXTURES ~ low-level               }....................
#FIXME: This fixture should no longer be required, suggesting that this entire
#module should be safely excisable. Consider it, please.
@fixture(scope='session')
def betse_init() -> None:
    '''
    Fixture initializing BETSE.

    This fixture permits BETSE functionality (e.g., classes, functions) to be
    safely imported and called by more fine-grained fixtures and tests. Failing
    to require this fixture _will_ reliably raise exceptions on importing or
    calling BETSE functionality.
    '''

    # Defer heavyweight imports to their point of use.
    from betse import ignition

    # Initialize BETSE.
    ignition.init()
