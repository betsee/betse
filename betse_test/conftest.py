#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
**Global test configuration** (i.e., early-time configuration guaranteed to be
run by :mod:`pytest` *after* passed command-line arguments are parsed) for
this test suite.

:mod:`pytest` implicitly imports *all* functionality defined by this module
into *all* submodules of this subpackage.

See Also
----------
:mod:`conftest`
    Root test configuration applied before this configuration.
'''

# ....................{ IMPORTS                           }....................
import pytest

# ....................{ IMPORTS ~ fixture : manual        }....................
# Import fixtures required to be manually required by other fixtures and tests.

from betse.util.test.pytest.fixture.pytfixture import monkeypatch_session
from betse_test.fixture.tempdirer import betse_temp_dir
from betse_test.fixture.simconf.simconfer import (
    betse_sim_conf,
    betse_sim_conf_default,
    betse_sim_conf_compat,
)

# ....................{ IMPORTS ~ fixture : autouse       }....................
# Import fixtures automatically run at the start of the current test session,
# typically *NOT* manually required by specific tests, *AFTER* importing all
# non-autouse fixtures possibly required by these autouse fixtures above.

#FIXME: This fixture has been temporarily superceded by the
#betse_test.conftest.pytest_configure() and pytest_unconfigure() hooks until
#inevitably required yet again.
# from betse_test.fixture.initter import betse_init_session

# ....................{ HOOKS ~ configure                 }....................
def pytest_configure(config) -> None:
    '''
    Hook run immediately *after* parsing all command-line options and loading
    all third-party :mod:`pytest` plugins (including application-specific
    ``conftest`` scripts) but *before* performing test collection.

    Specifically, this hook (in no particular order):

    * Instantiates and initializes the application metadata singleton in a
      manner suitable for unit testing.

    See Also
    ----------
    :func:`_init_app`
        Further details on application initialization.
    '''

    # Defer heavyweight imports.
    from betse_test.fixture import initter

    # Prepend a leading newline, which py.test curiously neglects to do itself.
    print('\n')

    # Initialize the application metadata singleton.
    initter.init_app()


def pytest_unconfigure(config) -> None:
    '''
    Hook run immediately *before* exiting the current test session.
    '''

    # Defer heavyweight imports.
    from betse.util.test.pytest import pytests
    from betse_test.fixture import initter

    # Deinitialize the application metadata singleton.
    initter.deinit_app()

    # Print this closure.
    pytests.output('Ego, ergo simulare.')

# ....................{ HOOKS ~ test                      }....................
def pytest_runtest_setup(item: 'pytest.main.Item') -> None:
    '''
    Hook run immediately *before* running the passed test.

    Specifically:

    * If this is a **serial test** (i.e., test method bound to an instance of
      the :class:`betse.util.test.pytest.pytabc.SerialTestABC` superclass) for
      which a prior serial test in the same test class was recorded as failing
      by the :func:`pytest_runtest_makereport` hook, this hook marks this test
      as xfailing.

    Parameters
    ----------
    item : pytest.main.Item
        Metadata encapsulating this test callable (e.g., function, method).

    See Also
    ----------
    https://pytest.org/latest/example/simple.html#incremental-testing-test-steps
        Official py.test code snippet inspiring this implementation.
    '''

    # Defer heavyweight imports.
    from betse.exceptions import BetseTestHookException
    from betse.util.test.pytest.pytabc import SerialTestABC

    # For each list of fixtures requested by this test...
    for fixture_defs in item._fixtureinfo.name2fixturedefs.values():
        # For each such fixture...
        for fixture_def in fixture_defs:
            # True only if this fixture is parametrized serially (i.e., is
            # decorated by the @serialize_parametrized_fixture decorator).
            is_fixture_parametrized_serially = getattr(
                fixture_def.func,
                '_betse_is_fixture_parametrized_serially',
                False)
            # print('funcinfo: {}'.format(fixture_def))
            # print('funcinfo.func: {}'.format(fixture_def.func))
            # print('funcinfo.func.is: {}'.format(is_fixture_parametrized_serially))

            # If this fixture is parametrized serially, propagate this setting
            # to this test by dynamically marking this test accordingly.
            if is_fixture_parametrized_serially:
                item.keywords['serialize_parametrized_test'] = True

    # If this test is parametrized serially, do so *AFTER* propagating this
    # setting from fixtures requested by this test above.
    if 'serialize_parametrized_test' in item.keywords:
        # If this test is *NOT* parametrized, raise an exception.
        if not hasattr(item, '_genid'):
            raise BetseTestHookException(
                'Unparametrized test "{}" marked as serial.'.format(item.name))

        # Test callable (i.e., underlying function or method object).
        test_callable = item.obj

        # Unique identifier of the first failing parameters passed to this test
        # if any or "None" otherwise.
        first_failing_param_id = getattr(
            test_callable, '_betse_first_failing_param_id', None)

        # If any parameters passed to this test failed, mark these subsequently
        # passed parameters as also xfailing.
        if first_failing_param_id is not None:
            pytest.xfail(
                'Prior serial test parametrization "{}" failed.'.format(
                    first_failing_param_id))
    #FIXME: This should *REALLY* simply be implemented as a new
    #@serialize_class decorator, much like @serialize_parametrized_test above.

    # If this is an unparametrized test intended to be run serially, do so.
    elif SerialTestABC.is_test_serial(item):
        # Object to which this test method is bound.
        test_instance = item.parent

        # Name of the first previously declared test method to fail in this
        # test method's class if any or "None" otherwise.
        first_failure_method_name = getattr(
            test_instance, '_first_failure_method_name', None)

        # If a previously declared test failed, mark this subsequently declared
        # test as also xfailing.
        if first_failure_method_name is not None:
            pytest.xfail('Prior serial test "{}" failed.'.format(
                test_instance._first_failure_method_name))


def pytest_runtest_makereport(
    item: 'pytest.main.Item', call: 'pytest.runner.CallInfo'):
    '''
    Hook run immediately *after* the passed test returned the passed result.

    Specifically:

    * If this is a **serial test** (i.e., test method bound to an instance of
      the the :class:`betse.util.test.pytest.pytabc.SerialTestABC` superclass)
      that failed, this hook records this failure for subsequent analysis by
      the :func:`pytest_runtest_setup` hook.

    Parameters
    ----------
    item : pytest.main.Item
        Metadata encapsulating this test callable (e.g., function, method).
    call : pytest.runner.CallInfo
        Metadata encapsulating the value returned or exception raised by
        calling this test callable.

    See Also
    ----------
    https://pytest.org/latest/example/simple.html#incremental-testing-test-steps
        Official py.test code snippet inspiring this implementation.
    '''

    # Defer heavyweight imports.
    from betse.util.test.pytest.pytabc import SerialTestABC

    # If this test failed...
    if call.excinfo is not None:
        # If this is a parametrized test intended to be run serially...
        if 'serialize_parametrized_test' in item.keywords:
            # Test callable (i.e., underlying function or method object).
            test_callable = item.obj

            # Record the unique identifiers of these test parameters for
            # subsequent use by the pytest_runtest_setup() hook.
            test_callable._betse_first_failing_param_id = item._genid
        #FIXME: This should *REALLY* simply be implemented as a new
        #@class_serial decorator, much like @parametrize_serial above.
        # If this is an unparametrized test intended to be run serially...
        elif SerialTestABC.is_test_serial(item):
            # Object to which this test method is bound.
            test_parent = item.parent

            # Record the name of this test for subsequent use by the
            # pytest_runtest_setup() hook.
            test_parent._first_failure_method_name = item.name
