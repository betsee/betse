#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Decorators marking tests and fixtures as being **parametrized** (i.e., accepting
two or more sets of parameters passed to these tests and fixtures).

Most of these decorators mark tests and fixtures with BETSE-specific keywords
inspected _only_ by BETSE-specific py.test hooks defined by `conftest` plugins.
'''

# ....................{ IMPORTS                            }....................
import pytest
from betse.util.type.types import type_check, CallableTypes
from collections import OrderedDict
# from functools import wraps

# ....................{ PARAMS                             }....................
# Outer decorator accepting all parameters explicitly passed to this decorator.
# Sadly, these parameters cannot simply be passed as an ordered dictionary of
# variadic keyword arguments (i.e., order-preserving "**kwargs"). Since PEP 468
# has yet to be accepted, Python has yet to support such arguments. See:
#     http://legacy.python.org/dev/peps/pep-0468
@type_check
def parametrize_test(param_name_to_values: OrderedDict) -> CallableTypes:
    '''
    Parametrize the decorated test callable with the passed ordered dictionary
    mapping the name to value of each parameter to pass to that test (in order).

    Fixtures
    ----------
    This decorator should _not_ be applied to fixtures, for which the
    `@pytest.fixture()` decorator should be applied instead as follows:

    * Pass that decorator:
    * The mandatory `params` keyword argument, whose value is the list of all
      sets of parameters to be iteratively passed to that fixture.
    * The optional `ids` keyword argument, whose value is the list of all
      human-readable unique identifiers to be assigned to each such set of
      parameters (in the same order). While optional, this argument premains
      highly recommended.
    * Pass that fixture the `request` builtin fixture, whose `param` attribute
      supplies the set of parameters passed to the current fixture invocation.

    Parameters
    ----------
    param_name_to_values : OrderedDict
        Ordered dictionary whose:
        * Keys are the Python-conformant names of all parameters to be passed to
          the decorated test.
        * Values are sequences (e.g., lists, tuples) of the values of those
          parameters to be iteratively passed, where the:
          . First element of each such sequence is the first value of that
            parameter to be passed on the first parametrization of that test.
          . And so on.
    '''

    # Inner closure decorating the actual test callable.
    def _parametrize_test_inner(test_callable) -> CallableTypes:
        # Defer to the canonical parametrize() decorator.
        return pytest.mark.parametrize(
            # Comma-delimited string listing the names of these parameters,
            # converted from this dictionary's keys.
            argnames=', '.join(param_name_to_values.keys()),

            # Sequence of the values of each parametrization of these
            # parameters.
            argvalues=param_name_to_values.values(),
        )(test_callable)

    # Return this inner closure.
    return _parametrize_test_inner

# ....................{ PARAMS ~ serial                    }....................
serialize_parametrized_test = pytest.mark.serialize_parametrized_test
'''
Mark the decorated parametrized test as **serial** (i.e., parametrized such that
the success of each subsequent parameter set depends on the success of all
previous parameter sets for the current test requesting this decoration).

If this test is _not_ parametrized, an exception is raised. Else, each parameter
set passed to this test is tested serially. On the first failure of this test
passed a parameter set, each subsequent call to this test passed a subsequent
parameter set will be automatically marked as an `XFAIL` and hence fail
_without_ being run.

Implementation
----------
The majority of the black magic required by this decoration is implemented as
low-level `py.test` hooks in the top-level :mod:`betse_test.conftest` plugin. To
preserve state between parametrized calls to the same test, these hooks
dynamically add the following BETSE-specific attributes to this test's
underlying function or method object:

* `_betse_first_failing_param_id`, the unique identifier of the first set of
  parameter values passed to this test raising an exception for the current test
  session if any _or_ `None` otherwise (i.e., if this test has yet to raise an
  exception for any parameters).
'''


def serialize_parametrized_fixture(fixture):
    '''
    Decorator marking the decorated parametrized fixture as **serial** (i.e.,
    parametrized such that the success of each subsequent parameter set depends
    on the success of all previous parameter sets for the current test
    requesting this decoration).

    If this fixture is _not_ parametrized, an exception is raised.

    This decoration is automatically propagated to all tests requesting this
    fixture. Decorating a fixture with this decorator is functionally
    equivalent to decorating all tests requesting this fixture by the
    `@serialize_parametrized_test` decorator.

    Implementation
    ----------
    Technically, `py.test` currently ignores marks decorating fixtures. To
    circumvent this omission, this decorator adds the following BETSE-specific
    attribute to this fixture's underlying function or method object:

    * `_betse_is_fixture_parametrized_serially`, `True` only if this fixture is
      parametrized serially. If this fixture is _not_ parametrized serially,
      this attribute is usually undefined.

    See Also
    ----------
    :func:`serialize_parametrized_test`
        Further details on serial parametrization.
    '''

    # Mark this fixture as serially parametrized. Since fixtures do *NOT*
    # support marks, a BETSE-specific attribute is added to this fixturegg
    fixture._betse_is_fixture_parametrized_serially = True
    return fixture
