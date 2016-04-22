#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Placeholder functional tests.
'''

# ....................{ IMPORTS                            }....................
# from betse.util.py import callers
# from betse.util.io import stdouts
# from test.func.conftest.sim_context import SimTestContext
# from func.conftest.sim_context import SimTestContext
from betse_test_func.util.context import fixture
import pytest

# ....................{ FIXTURES                           }....................
@pytest.fixture()
def _sub_fixture(request):
    # print('caller: {}'.format(callers.get_caller_basename()))
    # stdouts.output_traceback()
    print('fixturenames: {}'.format(request.fixturenames))

@pytest.fixture()
def null_fixture(_sub_fixture):
    pass

# ....................{ COMMANDS                           }....................
def test_cli_failure(tmpdir, null_fixture) -> None:
    '''
    Placeholder failing test.
    '''
    raise ValueError('Placeholder failing test.')

def test_cli_success(null_fixture) -> None:
    '''
    Placeholder successful test.
    '''
    pass
