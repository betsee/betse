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

    # If this test callable is a method intended to be run serially (i.e., in
    # test method declaration order such that subsequently declared test methods
    # depend on the success of all previously declared test methods)...
    if SerialTestABC.is_test_serial(item):
        # Object to which this test method is bound.
        test_instance = item.parent

        # If a previously declared test method in this test method's class
        # failed, mark this subsequently declared test as xfailing.
        if test_instance._first_failure_method_name is not None:
            pytest.xfail('prior serial test failed ({})'.format(
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

    # If this test callable is a method intended to be run serially (i.e., in
    # test method declaration order such that subsequently declared test methods
    # depend on the success of all previously declared test methods)...
    if SerialTestABC.is_test_serial(item):
        # If this test failed...
        if call.excinfo is not None:
            # Object to which this test method is bound.
            test_instance = item.parent

            # Record the name of this test method for subsequent use by the
            # pytest_runtest_setup() hook.
            test_instance._first_failure_method_name = item.name
