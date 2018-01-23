#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Decorators marking tests and fixtures as being **parametrized** (i.e., accepting
two or more sets of parameters passed to these tests and fixtures).

Most of these decorators mark tests and fixtures with BETSE-specific keywords
inspected *only* by BETSE-specific py.test hooks defined by ``conftest``
plugins.
'''

# ....................{ IMPORTS                            }....................
import itertools, pytest
from betse_test.exceptions import BetseTestParamException
from betse.util.type.types import (
    type_check, CallableTypes, ContainerType, MappingType, SequenceTypes,)

# ....................{ PARAMS                             }....................
# Outer decorator accepting all parameters explicitly passed to this decorator.
@type_check
def parametrize_test_setwise(
    params: MappingType = None,
    fixtures: MappingType = None,
    ids: SequenceTypes = None,
) -> CallableTypes:
    '''
    Optionally directly parametrize all non-fixture parameters accepted by the
    decorated test callable with the first passed dictionary _and_ optionally
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
    params : optional[MappingType]
        Dictionary mapping from the names of non-fixture parameters accepted by
        this test to the sequences of all values of those parameters to be
        iteratively passed to this test (_in order_) such that:
        . The first element of each such sequence is the first value of the
          corresponding parameter to be passed to the first parametrization of
          this test.
        . The second element of the same sequence is the second value of this
          parameter to be passed to the second parametrization of this test.
        . And so on.
        Defaults to `None`, in which case this test is assumed to accept no
        non-fixture parameters.
    fixtures : optional[MappingType]
        Dictionary mapping from the names of fixture parameters accepted by this
        test to the sequences of all values of those parameters to be
        iteratively passed to this test (_in order_) such that the same
        interpretation as for `params` holds. Defaults to `None`, in which case
        this test is assumed to accept no fixture parameters.
    ids : optional[SequenceTypes]
        Sequence of all human-readable labels uniquely identifying each
        parametrization of this test (_in the same order_). Defaults to `None`,
        in which non-human-readable labels will be automatically synthesized
        from the values of the parameters comprising these parametrizations.

    Raises
    ----------
    BetseTestParamException
        If the `ids` argument identifies a different number of parametrizations
        than the `params` and `fixtures` arguments actually passed.

    Examples
    ----------
    >>> from pytest import fixture
    >>> from betse_test.util.mark.param import parametrize_test_setwise
    >>> @fixture
    ... def feathered_serpent(serpent_name: str) -> str:
    ...     if serpent_name == 'Hualpa'
    ...         return 'Amazonia'
    ...     elif serpent_name == 'Mujaji'
    ...         return 'Azania'
    >>> @parametrize_test_setwise(
    ...     params={
    ...         'western_dragon': ('Celedyr', 'Hestaby',),
    ...         'eastern_dragon': ('Masaru', 'Ryumyo',),
    ...     },
    ...     fixtures={
    ...         'feathered_serpent': ('Hualpa', 'Mujaji',),
    ...     },
    ...     ids=('bad-dragons', 'good-dragons',),
    ... )
    ... def test_great_dragons(
    ...     western_dragon: str, eastern_dragon: str, feathered_serpent: str):
    ...     assert western_dragon in ('Celedyr', 'Hestaby',)
    ...     assert eastern_dragon in ('Masaru', 'Ryumyo',)
    ...     assert feathered_serpent in ('Amazonia', 'Azania',)
    '''

    # Defer heavyweight imports.
    from betse.util.type.iterables import zip_isometric

    # Default unpassed arguments to the empty dictionary for sanity. For unknown
    # (but presumably valid) reasons, these defaults *MUST* be established here
    # rather than in the inner closure defined below. Just. Because.
    if params   is None: params   = {}
    if fixtures is None: fixtures = {}

    # Inner closure decorating the actual test callable.
    def _parametrize_test_fixtures_inner(test_callable) -> CallableTypes:
        # Validate all non-fixture parametrizations to be sequences.
        for param_values in params.values():
            if not isinstance(param_values, SequenceTypes):
                raise BetseTestParamException(
                    'Parametrized test "{}" '
                    'non-fixture parametrization {} not a sequence.'.format(
                    test_callable.__name__, param_values))

        # Validate all fixture parametrizations to be sequences.
        for fixture_values in fixtures.values():
            if not isinstance(fixture_values, SequenceTypes):
                raise BetseTestParamException(
                    'Parametrized test "{}" '
                    'fixture parametrization {} not a sequence.'.format(
                    test_callable.__name__, fixture_values))

        # Tuple of the names of all fixture parameters (in arbitrary order) if
        # any or False otherwise. Attempting to pass an empty tuple as the
        # "indirect" keyword argument to the @pytest.mark.parametrize
        # decorator raises the following non-human-readable exception:
        #
        #     /usr/lib64/python3.4/site-packages/_pytest/python.py:1952: in pytest_generate_tests
        #         argnames = func_params[0]
        #     E   IndexError: tuple index out of range
        #
        # Likewise, attempting to pass None results in even more obscure
        # non-human-readable exceptions resembling:
        #
        #     /usr/lib64/python3.4/site-packages/_pytest/python.py:869: in setmulti
        #         valtype_for_arg = valtypes[arg]
        #     E   KeyError: '{param_name}'
        #
        # ...where "{param_name}" is the name of the first parameter listed in
        # the "param_fixture_names" tuple defined below.
        #
        # py.test, you are nice; but your failure to raise human-readable errors
        # is a deplorable waste worthy of unclean roadhouse outhouses.
        fixture_names = tuple(fixtures.keys()) or False

        # Tuple of the names of both non-fixture and fixture parameters (in a
        # predictable, but arbitrary, order).
        #
        # Note that the dict.keys() method returns instances of the "dict_keys"
        # class; unlike standard sequence classes (e.g., "list", "tuple"), this
        # class does *NOT* overload the "+" operator. Aggregating multiple
        # instances of this class requires leveraging the generic
        # itertools.chain() function, accepting arbitrary iterables. (Ugh.)
        param_fixture_names = tuple(itertools.chain(
            params.keys(), fixtures.keys()))
        # param_fixture_names = ','.join(itertools.chain(
        #     params.keys(), fixtures.keys()))

        # Tuple of n-tuples of all values for each parametrization of these
        # non-fixture and parameters (in this same predictable order). See
        # parametrize_test_setwise() for further discussion.
        param_fixture_values = tuple(zip_isometric(*itertools.chain(
            params.values(), fixtures.values())))

        # If the caller identified a different number of parametrizations than
        # were actually passed, raise an exception.
        if ids is not None and len(ids) != len(param_fixture_values):
            raise BetseTestParamException(
                'Parametrized test "{}" '
                'number of parametrization identifiers {} differs from '
                'number of parametrizations {}.'.format(
                 test_callable.__name__, len(ids), len(param_fixture_values)))
        # print('\n    test: {}; param_names: {!r}; values: {!r}; fixtures: {!r}; ids: {!r}'.format(test_callable.__name__, param_fixture_names, param_fixture_values, fixture_names, ids))

        # Defer to the canonical parametrize() decorator.
        return pytest.mark.parametrize(
            argnames=param_fixture_names,
            argvalues=param_fixture_values,
            indirect=fixture_names,
            ids=ids,
        )(test_callable)

    # Return this inner closure.
    return _parametrize_test_fixtures_inner

# ....................{ PARAMS ~ setwise                   }....................
# Outer decorator accepting all parameters explicitly passed to this decorator.

#FIXME: Actually implement. See the "Examples" docstring section for the
#intended usage of this decorator; basically, it provides an alternative means
#(and arguably equally intuitive) of parametrizing tests to the above
#parametrize_test_setwise() decorator.
#FIXME: The principal disadvantage of both this and the above
#parametrize_test_setwise() decorator is that neither supports the fine-grained
#specification of specific parametrizations to be either skipped or xfailed, as
#supported by the lower-level @pytest.mark.parametrize decorator. See the
#following official py.test documentation for a useful example:
#    http://doc.pytest.org/en/latest/skipping.html#skip-xfail-with-parametrize
#While the above parametrize_test_setwise() decorator *CANNOT* (by design) be
#refactored to support such fine-grained specification, this decorator should be
#able to do so. How? By defining two new "ParametrizationSkipIf" and
#"ParametrizationXFail" classes instantiated as follows:
#
#    @parametrize_test_setwise(
#        parametrizations={
#            'bad-dragons': {
#                'western_dragon':    'Celedyr',
#                'eastern_dragon':    'Masaru',
#                'feathered_serpent': 'Hualpa',
#            },
#            'good-dragons': ParametrizationXFail(
#                parametrizations={
#                    'western_dragon':    'Hestaby',
#                    'eastern_dragon':    'Ryumyo',
#                    'feathered_serpent': 'Mujaji',
#                },
#                reason="Good dragons raise unexpected breath attack.",
#            },
#        },
#        fixture_names={'feathered_serpent',},
#    )
#    def test_great_dragons(
#        western_dragon: str, eastern_dragon: str, feathered_serpent: str):
#        assert western_dragon in ('Celedyr', 'Hestaby',)
#        assert eastern_dragon in ('Masaru', 'Ryumyo',)
#        assert feathered_serpent in ('Amazonia', 'Azania',)
#
#This decorator would then need to explicitly differentiate such objects from
#standard dictionaries and behave as expected. Quite a bit of work, sadly, but
#certainly feasible.
#FIXME: Repair documentation.

@type_check
def parametrize_test_paramwise(
    fixture_names: ContainerType = None,
    parametrizations: MappingType = None,
) -> CallableTypes:
    '''
    Optionally directly parametrize all non-fixture parameters accepted by the
    decorated test callable with the first passed dictionary _and_ optionally
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
    params : optional[MappingType]
        Dictionary mapping from the names of non-fixture parameters accepted by
        this test to the sequences of all values of those parameters to be
        iteratively passed to this test (_in order_) such that:
        . The first element of each such sequence is the first value of the
          corresponding parameter to be passed to the first parametrization of
          this test.
        . The second element of the same sequence is the second value of this
          parameter to be passed to the second parametrization of this test.
        . And so on.
        Defaults to `None`, in which case this test is assumed to accept no
        non-fixture parameters.
    fixtures : optional[MappingType]
        Dictionary mapping from the names of fixture parameters accepted by this
        test to the sequences of all values of those parameters to be
        iteratively passed to this test (_in order_) such that the same
        interpretation as for `params` holds. Defaults to `None`, in which case
        this test is assumed to accept no fixture parameters.
    ids : optional[SequenceTypes]
        Sequence of all human-readable labels uniquely identifying each
        parametrization of this test (_in the same order_). Defaults to `None`,
        in which non-human-readable labels will be automatically synthesized
        from the values of the parameters comprising each such parametrization.

    Raises
    ----------
    BetseTestParamException
        If the `ids` argument identifies a different number of parametrizations
        than the `params` and `fixtures` arguments actually passed.

    Examples
    ----------
    >>> from pytest import fixture
    >>> from betse_test.util.mark.param import parametrize_test_setwise
    >>> @fixture
    ... def feathered_serpent(serpent_name: str) -> str:
    ...     if serpent_name == 'Hualpa'
    ...         return 'Amazonia'
    ...     elif serpent_name == 'Mujaji'
    ...         return 'Azania'
    >>> @parametrize_test_setwise(
    ...     parametrizations={
    ...         'bad-dragons': {
    ...             'western_dragon':    'Celedyr',
    ...             'eastern_dragon':    'Masaru',
    ...             'feathered_serpent': 'Hualpa',
    ...         },
    ...         'good-dragons': {
    ...             'western_dragon':    'Hestaby',
    ...             'eastern_dragon':    'Ryumyo',
    ...             'feathered_serpent': 'Mujaji',
    ...         },
    ...     },
    ...     fixture_names={'feathered_serpent',},
    ... )
    ... def test_great_dragons(
    ...     western_dragon: str, eastern_dragon: str, feathered_serpent: str):
    ...     assert western_dragon in ('Celedyr', 'Hestaby',)
    ...     assert eastern_dragon in ('Masaru', 'Ryumyo',)
    ...     assert feathered_serpent in ('Amazonia', 'Azania',)
    '''

    raise ValueError('This decorator currently unimplemented.')

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
