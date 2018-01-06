#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Fixtures creating temporary directories isolated for safety to specific tests.
'''

# ....................{ IMPORTS                            }....................
from betse_test.util import requests
from pytest import fixture

# ....................{ FIXTURES                           }....................
# Test-scope fixture creating and returning a new object for each discrete test.
@fixture
def betse_temp_dir(
    request: '_pytest.python.FixtureRequest',
    tmpdir_factory: '_pytest.tmpdir.tmpdir_factory',
) -> 'LocalPath':
    '''
    Per-test fixture creating a temporary directory isolated to the test
    requesting this fixture and returning an object encapsulating this
    directory.

    This directory is guaranteed to be specific to this test. The basename of
    this directory is the name of this test excluding the prefixing substring
    ``test_``. When requested by the ``test_cli_sim_default`` test, for example,
    this fixture creates a temporary directory ``{tmpdir}/cli_sim_default`` for
    the absolute path ``{tmpdir}`` of this test session's root temporary
    directory (e.g., ``/tmp/pytest-0/cli_sim_default``).

    This directory is safely accessible *only* for the duration of this test.
    Subsequently run tests and fixtures *cannot* safely reuse this directory,
    although doing so is technically feasible in unreliable edge-cases.

    Parameters
    ----------
    request : _pytest.python.FixtureRequest
        Builtin fixture describing the parent fixture or test of this fixture.
    tmpdir_factory : _pytest.tmpdir.tmpdir_factory
        Builtin session-scoped fixture whose ``mktemp()`` method returns a
        :class:`py.path.local` instance encapsulating a new temporary directory.

    Returns
    ----------
    LocalPath
        Object encapsulating this temporary directory.
    '''

    # Defer heavyweight imports.
    from betse.util.type.text import strs

    # Name of the current test.
    test_name = requests.get_tested_name(request)
    # print('    request.node: {}'.format(request.node))
    # print('    test_name: {}'.format(test_name))

    # Basename of the temporary directory containing this configuration file,
    # set to the name of the current test excluding the prefixing "test_".
    temp_dir_basename = strs.remove_prefix(
        text=test_name,
        prefix='test_',
        exception_message=(
            'Test name "{}" not prefixed by "test_".'.format(test_name)),
    )

    # Create this temporary directory and wrap this directory's absolute path
    # with a high-level "py.path.local" object. See also:
    #     http://pytest.org/latest/tmpdir.html#the-tmpdir-factory-fixture
    temp_dirpath = tmpdir_factory.mktemp(temp_dir_basename)

    # Return this object.
    return temp_dirpath
