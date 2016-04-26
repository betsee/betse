#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Unit tests for BETSE's Python API testing the importability of fragile public
modules and packages (e.g., containing cyclic import statements subject to
breaking by brute-force PyCharm refactorings).
'''

# ....................{ IMPORTS                            }....................
from betse_test.mark.skip import skip_unless_plugin_xdist

# ....................{ TESTS                              }....................
@skip_unless_plugin_xdist
def test_import_logs() -> None:
    '''
    Test the importability of BETSE's logging API.

    This API contains cyclic import statements commonly broken by PyCharm's
    brute-force refactorings.
    '''

    # Import this API and call all functionality containing these statements.
    from betse.util.io.log import logs
    logs.get()
