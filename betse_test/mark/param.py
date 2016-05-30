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
# from betse.util.type import types

# ....................{ PARAMS ~ alias                     }....................
parametrize_test = pytest.mark.parametrize
'''
Parametrize the decorated test callable, class, or module with the passed lists
of parameter names and values.

Fixtures
----------
This decorator should _not_ be applied to fixtures, for which the
`@pytest.fixture()` decorator should be applied instead as follows:

* Pass that decorator:
  * The mandatory `params` keyword argument, whose value is the list of all sets
    of parameters to be iteratively passed to that fixture.
  * The optional `ids` keyword argument, whose value is the list of all
    human-readable unique identifiers to be assigned to each such set of
    parameters (in the same order). While optional, this is highly recommended.
* Pass that fixture the `request` builtin fixture, whose `param` attribute
  supplies the set of parameters passed to the current fixture invocation.

Parameters
----------
argnames : str
    Comma-delimited string listing the Python-conformant names of all parameters
    to be iteratively passed to the decorated test.
argvalues : list
    List of all sets of parameters to be iteratively passed to the decorated
    test.
'''

# ....................{ PARAMS ~ serial                    }....................
parametrize_test_serial = pytest.mark.parametrize_test_serial
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
low-level py.test hooks in the top-level `betse_func.conftest` plugin. To
preserve state between parametrized calls to the same test, these hooks
dynamically add the following BETSE-specific attributes to this test's
underlying function or method object:

* `_betse_first_failing_param_id`, the unique identifier of the first set of
  parameter values passed to this test raising an exception for the current test
  session if any _or_ `None` otherwise (i.e., if this test has yet to raise an
  exception for any parameters).
'''


def parametrize_fixture_serial(fixture):
    '''
    Decorator marking the decorated parametrized fixture as **serial** (i.e.,
    parametrized such that the success of each subsequent parameter set depends
    on the success of all previous parameter sets for the current test
    requesting this decoration).

    If this fixture is _not_ parametrized, an exception is raised.

    This decoration is automatically propagated to all tests requesting this
    fixture. Decorating a fixture with this decorator is functionally
    equivalent to decorating all tests requesting this fixture by the
    `@parametrize_test_serial` decorator.

    Implementation
    ----------
    Technically, `py.test` currently ignores marks decorating fixtures. To
    circumvent this omission, this decorator adds the following BETSE-specific
    attribute to this fixture's underlying function or method object:

    * `_betse_is_fixture_parametrized_serial`, `True` only if this fixture is
      parametrized serially. If this fixture is _not_ parametrized serially,
      this attribute is usually undefined.

    See Also
    ----------
    `@parametrize_test_serial`
        Further details on serial parametrization.
    '''

    # Mark this fixture as serially parametrized. Since fixtures do *NOT*
    # support marks, a BETSE-specific attribute is added to this fixturegg
    fixture._betse_is_fixture_parametrized_serial = True
    return fixture
