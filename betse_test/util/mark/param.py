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
from betse.util.type.types import CallableTypes
from betse_test.exceptions import BetseTestParameterException
# from functools import wraps

# ....................{ PARAMS                             }....................
# Outer decorator accepting all parameters explicitly passed to this decorator.
# Sadly, these parameters cannot simply be passed as an ordered dictionary of
# variadic keyword arguments (i.e., order-preserving "**kwargs"). Since PEP 468
# has yet to be accepted, Python has yet to support such arguments. See:
#     http://legacy.python.org/dev/peps/pep-0468

#FIXME: Revise docstring.
def parametrize_test(**param_name_to_values) -> CallableTypes:
    '''
    Parametrize the decorated test callable with the passed tuple of parameter
    names and values to iteratively pass to that test (_in order_).

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
    param_name_values : tuple
        Tuple of sequential parameter names and keys to parametrize this test
        with. Each element of this tuple with:
        * Even index (e.g., the first and third elements) specifies the name of
          an existing parameter accepted by this test to be parametrized.
        * Odd index (e.g., the second and fourth elements) specifies a sequence
          of one or more values for the parameter whose name is specified by the
          preceding element, where the:
          . First element of this sequence is the first value of this parameter
            to be passed to the first parametrization of this test.
          . And so on.

    Raises
    ----------
    :exc:`betse_test.exceptions.BetseTestParameterException`
            If this tuple is _not_ of even length.
    '''

    # Defer heavyweight imports.
    from betse.util.type import ints

    # Inner closure decorating the actual test callable.
    def _parametrize_test_inner(test_callable) -> CallableTypes:
        # Comma-delimited string listing the names of these parameters.
        param_names = ','.join(param_name_to_values.keys())

        # Tuple of n-tuples of all values for each parametrization of these
        # parameters. Curiously, py.test has been hard-coded to require a
        # statically sized container rather than a dynamically sized generator.
        # Failure to reduce this zip() object to a tuple induces py.test to
        # ignore all tests decorated with this decorator with a warning
        # resembling:
        #
        #    SKIP [1] /usr/lib64/python3.4/site-packages/_pytest/python.py:1417:
        #    got empty parameter set, function test_params at
        #    /home/diogenes/py/betse/betse_test/test_params.py:13
        #
        # Dismantled, this logic is:
        #
        # * "param_name_values[1::2]", a tuple of each sequence of values
        #   specific to each parameter.
        # * "*", unpacking these sequences out of this tuple.
        # * zip(...), repacking the same element of each such sequence into a
        #   new sequence of all such elements encapsulated by this zip object.
        # * tuple(...), converting this zip object into a tuple.
        #
        # Note that the official Python documentation explicitly guarantees the
        # order of elements returned by the keys() and values() methods of
        # unmodified dictionaries to correspond:
        #
        #    "If items(), keys(), values(), iteritems(), iterkeys(), and
        #     itervalues() are called with no intervening modifications to the
        #     dictionary, the order of items will directly correspond."
        #
        # zip: obfuscating Python since 1989.
        param_values = tuple(zip(*param_name_to_values.values()))
        # print('\n!!!!!param_names: {!r}; values: {!r}'.format(param_names, list(param_values)))

        #FIXME: Raise an exception unless *ALL* parameter values sequences are
        #of the same length. To do so sanely, leverage the newly defined
        #zip_isomorphic() utility generator.

        # Defer to the canonical parametrize() decorator.
        return pytest.mark.parametrize(
            argnames=param_names, argvalues=param_values,
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
