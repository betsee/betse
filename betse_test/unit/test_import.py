#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Unit tests for BETSE's Python API testing the importability of fragile public
modules and packages (e.g., containing cyclic import statements subject to
breaking by brute-force PyCharm refactorings).
'''

# ....................{ IMPORTS                            }....................
from betse_test.util.mark.fail import xfail
# from betse_test.mark.skip import skip_unless_plugin_xdist

# ....................{ TESTS                              }....................
#FIXME: Ugh! Ideally, this test should *NOT* require xdist support. This test
#must be isolated to a subprocess. To do so, the simplest means is the xdist
#"--boxed" CLI option. However, we'd only want *THIS* test to be isolated in
#that manner -- not *EVERY* test! To do so, we'll probably want to examine the
#xdist codebase, ascertain how "--boxed" is implemented, and provide a new
#utility test function or decorator performing the equivalent. Alternately,
#PyInstaller's testing framework already supports test isolation (probably using
#the "subprocess" module), suggesting we could also examine that as a fallback.

# @skip_unless_plugin_xdist
@xfail(reason='Test isolation unsupported.')
def test_import_logs() -> None:
    '''
    Test the importability of BETSE's logging API.

    This API contains cyclic import statements commonly broken by PyCharm's
    brute-force refactorings.
    '''

    # Import this API and call all functionality containing these statements.
    from betse.util.io.log import logs
    logs.get()
