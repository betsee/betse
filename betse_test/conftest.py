#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Global test configuration for all tests.

`py.test` implicitly imports _all_ functionality defined by this module into
_all_ test modules.
'''

# ....................{ IMPORTS                            }....................
import pytest
from betse_test.util.exceptions import BetseTestHookException
from betse_test.util.testabc import SerialTestABC

# ....................{ IMPORTS ~ fixture                  }....................
from betse_test.fixture.ignition import betse_init

# ....................{ HOOKS                              }....................
def pytest_runtest_setup(item: 'pytest.main.Item'):
    '''
    Hook run immediately _before_ running the passed test.

    If this is a **serial test** (i.e., test method bound to an instance of the
    the `SerialTestABC` superclass) for which a prior serial test in the same
    test class was recorded as failing by the `pytest_runtest_makereport()`
    hook, this hook marks this test as xfailing.

    Parameters
    ----------
    item : pytest.main.Item
        Metadata encapsulating this test callable (e.g., function, method).

    See Also
    ----------
    https://pytest.org/latest/example/simple.html#incremental-testing-test-steps
        Official py.test code snippet inspiring this implementation.
    '''

    # If this is a parametrized test intended to be run serially, do so.
    if 'serial_parametrized' in item.keywords:
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
    #FIXME: This should *REALLY* simply be implemented as a new @serial
    #decorator, much like @serial_parametrized above.
    # If this is an unparametrized test intended to be run serially, do so.
    elif SerialTestABC.is_test_serial(item):
        # Object to which this test method is bound.
        test_instance = item.parent

        # Name of the first previously declared test method to fail in this test
        # method's class if any or "None" otherwise.
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
    Hook run immediately _after_ the passed test returned the passed result.

    If this is a **serial test** (i.e., test method bound to an instance of the
    the `SerialTestABC` superclass) that failed, this hook records this failure
    for subsequent analysis by the `pytest_runtest_setup()` hook.

    Parameters
    ----------
    item : pytest.main.Item
        Metadata encapsulating this test callable (e.g., function, method).
    call : pytest.runner.CallInfo
        Metadata encapsulating the value returned or exception raised by calling
        this test callable.

    See Also
    ----------
    https://pytest.org/latest/example/simple.html#incremental-testing-test-steps
        Official py.test code snippet inspiring this implementation.
    '''

    # If this test failed...
    if call.excinfo is not None:
        # If this is a parametrized test intended to be run serially...
        if 'serial_parametrized' in item.keywords:
            # Test callable (i.e., underlying function or method object).
            test_callable = item.obj

            # Record the unique identifiers of these test parameters for
            # subsequent use by the pytest_runtest_setup() hook.
            test_callable._betse_first_failing_param_id = item._genid
        #FIXME: This should *REALLY* simply be implemented as a new @serial
        #decorator, much like @serial_parametrized above.
        # If this is an unparametrized test intended to be run serially...
        elif SerialTestABC.is_test_serial(item):
            # Object to which this test method is bound.
            test_parent = item.parent

            # Record the name of this test for subsequent use by the
            # pytest_runtest_setup() hook.
            test_parent._first_failure_method_name = item.name
