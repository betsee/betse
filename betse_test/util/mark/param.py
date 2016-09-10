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
import itertools, pytest
from betse_test.exceptions import BetseTestParamException
from betse.util.type.types import (
    type_check, CallableTypes, MappingType, SequenceTypes,)

# ....................{ PARAMS                             }....................
# Outer decorator accepting all parameters explicitly passed to this decorator.
@type_check
def parametrize_test(
    params: MappingType,
    fixtures: MappingType = None,
    ids: SequenceTypes = None,
) -> CallableTypes:
    '''
    Directly parametrize all non-fixture parameters accepted by the decorated
    test callable with the first passed dictionary as well as optionally
    indirectly parametrize all fixture parameters accepted by this callable with
    the second passed dictionary.

    Fixtures
    ----------
    This decorator should _not_ be directly applied to fixtures, for which the
    :func:`pytest.fixture` decorator should be applied instead as follows:

    * Pass that decorator:
      * The mandatory `params` keyword argument, whose value is the sequence of
        all sequences of parameters to be iteratively passed to that fixture.
      * The optional `ids` keyword argument, whose value is the sequence of all
        human-readable unique identifiers to be assigned to each such set of
        parameters (in the same order). While optional, this argument premains
        highly recommended.
    * Pass that fixture the `request` builtin fixture, whose `param` attribute
      supplies the set of parameters passed to the current fixture invocation.

    Fixture Parameters
    ----------
    For each fixture accepted by the decorated test callable, the
    :func:`pytest.fixture` decorator should be applied to that fixture as
    follows:

    * Pass that fixture the `request` builtin fixture, whose `param` attribute
      supplies the set of parameters passed to the current fixture invocation.

    Avoid passing that decorator either the `params` or `ids` keyword argument,
    whose values will be indirectly supplied by this decorator to that fixture.

    Parameters
    ----------
    params : MappingType
        Dictionary mapping from the names of non-fixture parameters accepted by
        this test to the sequences of all values of those parameters to be
        iteratively passed to this test (_in order_) such that:
        . The first element of each such sequence is the first value of the
          corresponding parameter to be passed to the first parametrization of
          this test.
        . The second element of the same sequence is the second value of this
          parameter to be passed to the second parametrization of this test.
        . And so on.
    fixtures : optional[MappingType]
        Dictionary mapping from the names of fixture parameters accepted by this
        test to the sequences of all values of those parameters to be
        iteratively passed to this test (_in order_) such that the same
        interpretation as for `params` holds. Defaults to `None`, in which case
        this test is assumed to accept only non-fixture parameters.
    ids : optional[SequenceTypes]
        Sequence of all human-readable labels uniquely identifying each
        parametrization of this test (_in the same order_). Defaults to `None`,
        in which non-human-readable labels will be automatically synthesized
        from the values of the parameters comprising each such parametrization.

    Raises
    ----------
    :exc:`BetseTestParamException`
        If the `ids` argument identifies a different number of parametrizations
        than the `params` and `fixtures` arguments actually passed.

    Examples
    ----------
    >>> from pytest import fixture
    >>> from betse_test.util.mark.param import parametrize_test
    >>> @fixture
    ... def feathered_serpent(serpent_name: str) -> str:
    ...     if serpent_name == 'Hualpa'
    ...         return 'Amazonia'
    ...     elif serpent_name == 'Mujaji'
    ...         return 'Azania'
    >>> @parametrize_test(
    ...     params={
    ...         'western_dragon': ('Celedyr', 'Hestaby',),
    ...         'eastern_dragon': ('Masaru', 'Ryumyo',),
    ...     },
    ...     fixtures={
    ...         'feathered_serpent': ('Hualpa', 'Mujaji',),
    ...     },
    ...     ids=('bad-dragons', 'good-dragons',),
    ... )
    ... def test_fixture_params(
    ...     western_dragon: str, eastern_dragon: str, feathered_serpent: str):
    ...     assert western_dragon in ('Celedyr', 'Hestaby',)
    ...     assert eastern_dragon in ('Masaru', 'Ryumyo',)
    ...     assert feathered_serpent in ('Amazonia', 'Azania',)
    '''

    # Defer heavyweight imports.
    from betse.util.type.iterables import zip_isometric

    # If this test accepts no parametrized fixtures, default this argument to
    # the empty dictionary for sanity.
    if fixtures is None:
        fixtures = {}

    # Inner closure decorating the actual test callable.
    def _parametrize_test_fixtures_inner(test_callable) -> CallableTypes:
        # Tuple of the names of all fixture parameters (in arbitrary order).
        fixture_names = tuple(fixtures.keys())

        # Comma-delimited string listing the names of both non-fixture and
        # fixture parameters (in a predictable but arbitrary order).
        #
        # Note that the dict.keys() method returns instances of the "dict_keys"
        # class; unlike standard sequence classes (e.g., "list", "tuple"), this
        # class does *NOT* overload the "+" operator. Aggregating multiple
        # instances of this class requires leveraging the generic
        # itertools.chain() function, accepting arbitrary iterables. (Ugh.)
        param_names = ','.join(itertools.chain(params.keys(), fixtures.keys()))

        # Tuple of n-tuples of all values for each parametrization of these
        # non-fixture and parameters (in this same predictable order). See
        # parametrize_test() for further discussion.
        param_values = tuple(zip_isometric(*itertools.chain(
            params.values(), fixtures.values())))
        # print('\n!!!!!param_names: {!r}; values: {!r}; fixtures: {!r}'.format(param_names, param_values, fixture_names))

        # If the caller identified a different number of parametrizations than
        # were actually passed, raise an exception.
        if ids is not None and len(ids) != len(param_values):
            raise BetseTestParamException(
                'Number of parametrization identifiers {} differs from '
                'Number of parametrizations {}'.format(
                len(ids), len(param_values)))

        # Defer to the canonical parametrize() decorator.
        return pytest.mark.parametrize(
            argnames=param_names,
            argvalues=param_values,
            indirect=fixture_names,
            ids=ids,
        )(test_callable)

    # Return this inner closure.
    return _parametrize_test_fixtures_inner

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
