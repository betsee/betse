#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Utility functions simplifying access to :mod:`pytest` built-in ``request``
fixture.

The official documentation for this fixture is blatantly inadequate. Even were
this documentation adequate, however, this fixture's API is overly obtuse to an
absurd (and almost obfuscatory) degree. This module provides an intelligible
alternative to that API.

See Also
----------
https://pytest.org/latest/builtin.html#_pytest.python.FixtureRequest
    Official documentation for this fixture – such as it is.
'''

# Note that all functions defined below accepting a "request" fixture object
# assume the passed object to be of the expected type rather than validating
# this to be the case. Why?  Because this type resides in a private module that
# does *NOT* appear to be externally importable: e.g.,
#
#     >>> import _pytest.python.FixtureRequest as req
#     ImportError: cannot import name 'transfer_markers'

# ....................{ IMPORTS                            }....................
from pytest import Function
from betse.util.type import sequences
from betse.util.type.types import type_check, TestableTypes
from betse_test.exceptions import BetseTestFixtureException
from betse_test.util import pytests

# ....................{ CONSTANTS                          }....................
_DEFAULT_VALUE_NONE = object()
'''
Default value for the optional ``default_value`` argument accepted by the
:func:`get_fixture_param` function.

This value permits that function to distinguish between default values of value
``None`` and the lack of any default value altogether.
'''

# ....................{ EXCEPTIONS                         }....................
def die_unless_tested(request: '_pytest.python.FixtureRequest') -> None:
    '''
    Raise an exception unless the passed ``request`` fixture object was
    transitively requested by a test.

    Parameters
    ----------
    request : _pytest.python.FixtureRequest
        Object passed to fixtures and tests requesting the ``request`` fixture.

    Raises
    ----------
    BetseTestFixtureException
        If no test transitively requested this ``request`` fixture.

    See Also
    ----------
    :func:`is_tested`
        Further details.
    '''

    if not is_tested(request):
        raise BetseTestFixtureException(
            '"request" fixture transitively requested by no test '
            '(e.g., due to being requested by a file- or '
            'session-scoped fixture): {}'.format(request))

# ....................{ EXCEPTIONS ~ fixture               }....................
def die_unless_fixture(request: '_pytest.python.FixtureRequest') -> None:
    '''
    Raise an exception unless the passed ``request`` fixture object was directly
    requested by a fixture rather than a test.

    Parameters
    ----------
    request : _pytest.python.FixtureRequest
        Object passed to fixtures and tests requesting the ``request`` fixture.

    Raises
    ----------
    BetseTestFixtureException
        If a test rather than fixture requested this ``request`` fixture.
    '''

    if not is_fixture(request):
        raise BetseTestFixtureException(
            '"request" fixture requested by '
            'test rather than fixture: {}'.format(request))


def die_unless_fixture_parametrized(
    request: '_pytest.python.FixtureRequest') -> bool:
    '''
    Raise an exception unless the fixture requesting the passed ``request``
    fixture object was parametrized (either directly or indirectly).

    Parameters
    ----------
    request : _pytest.python.FixtureRequest
        Object passed to fixtures and tests requesting the ``request`` fixture.

    Raises
    ----------
    BetseTestFixtureException
        If either:
        * A test rather than fixture requested this ``request`` fixture.
        * A fixture that was _not_ parametrized requested this ``request``
          fixture.
    '''

    if not is_fixture_parametrized(request):
        raise BetseTestFixtureException(
            '"request" fixture requested by '
            'unparametrized fixture "{}".'.format(get_fixture_name(request)))

# ....................{ TESTERS                            }....................
def is_tested(request: '_pytest.python.FixtureRequest') -> bool:
    '''
    `True` only if the passed ``request`` fixture object was transitively
    requested by a test (e.g., directly by a test or indirectly by a test-scope
    fixture).

    Parameters
    ----------
    request : _pytest.python.FixtureRequest
        Object passed to fixtures and tests requesting the ``request`` fixture.

    Returns
    ----------
    bool
        `True` only if a test transitively requested this ``request`` fixture.
    '''

    return isinstance(request.node, Function)

# ....................{ TESTERS ~ fixture                  }....................
def is_fixture(request: '_pytest.python.FixtureRequest') -> bool:
    '''
    `True` only if the passed ``request`` fixture object was directly requested by
    a fixture rather than a test.

    Parameters
    ----------
    request : _pytest.python.FixtureRequest
        Object passed to fixtures and tests requesting the ``request`` fixture.

    Returns
    ----------
    bool
        `True` only if a fixture requested this ``request`` fixture.
    '''

    return request.fixturename is not None


def is_fixture_parametrized(request: '_pytest.python.FixtureRequest') -> bool:
    '''
    `True` only if the fixture requesting the passed ``request`` fixture object
    was parametrized (either directly or indirectly).

    Parameters
    ----------
    request : _pytest.python.FixtureRequest
        Object passed to fixtures and tests requesting the ``request`` fixture.

    Returns
    ----------
    bool
        `True` only if this fixture was parametrized.

    Raises
    ----------
    BetseTestFixtureException
        If a test rather than fixture requested this ``request`` fixture.
    '''

    # If this "request" fixture was requested by a test, raise an exception.
    die_unless_fixture(request)

    # Else, that fixture was itself requested by a fixture. Return whether on
    # not the latter was parametrized. py.test API: you fail again.
    return hasattr(request, 'param')

# ....................{ GETTERS                            }....................
def get_tested_name(request: '_pytest.python.FixtureRequest') -> str:
    '''
    Unqualified name of the current test (e.g., `test_cli_info`) if the passed
    ``request`` fixture object was either directly requested by that test *or*
    indirectly requested by a test-scope fixture.

    This function returns either:

    * If the current test is a test function, this function's unqualified name.
    * If the current test is a test method, this method's unqualified name.

    Parameters
    ----------
    request : _pytest.python.FixtureRequest
        Object passed to fixtures and tests requesting the ``request`` fixture.

    Returns
    ----------
    str
        Unqualified name of this test.

    Raises
    ----------
    BetseTestFixtureException
        If the passed ``request`` fixture object was *not* transitively requested
        by a test (e.g., was requested by a file- or session-scope fixture), in
        which case no test exists to obtain the name of.
    '''

    # Raise an exception unless a test transitively requested this fixture.
    die_unless_tested(request)

    # No one would ever think of trying this. Except one man did.
    return request.node.name

# ....................{ GETTERS ~ fixture                  }....................
def get_fixture_name(request: '_pytest.python.FixtureRequest') -> str:
    '''
    Name of the current fixture that explicitly requested the passed ``request``
    fixture object.

    If this ``request`` fixture was requested by a test rather than a fixture,
    an exception is raised.

    Parameters
    ----------
    request : _pytest.python.FixtureRequest
        Object passed to fixtures and tests requesting the ``request`` fixture.

    Returns
    ----------
    str
        Name of such fixture.
    '''

    # If this "request" fixture was requested by a test, raise an exception.
    die_unless_fixture(request)

    # Else, that fixture was itself requested by a fixture. Return the name of
    # the latter.
    return request.fixturename


@type_check
def get_fixture_param(
    request,
    default_value: object = _DEFAULT_VALUE_NONE,
    type_expected: TestableTypes = None,
) -> object:
    '''
    Parameter parametrizing the fixture requesting the passed ``request`` fixture,
    optionally defaulted to the passed default value if the former fixture is
    unparametrized _and_ optinionally validated to be an instance of the passed
    type or tuple of types.

    Parameters
    ----------
    request : _pytest.python.FixtureRequest
        Object passed to fixtures and tests requesting the ``request`` fixture.
    default_value: optional[object]
        Default value to return if this fixture is unparametrized. Defaults
        to :data:`_DEFAULT_VALUE_NONE`, in which case an exception is raised
        instead if this fixture is unparametrized.
    type_expected : optional[TestableTypes]
        Class or tuple of classes to validate this object to be an instance of.
        If this object is _not_ an instance of this class or classes, an
        exception is raised. Defaults to `None`, in which case this object is
        returned unvalidated.

    Returns
    ----------
    object
        Parameter parametrizing this fixture.

    Raises
    ----------
    BetseTestFixtureException
        If either:
        * A test rather than fixture requested this ``request`` fixture.
        * An unparametrized fixture requested this ``request`` fixture.
        * A fixture parametrized by an instance of a class _not_ the passed
          `check_type` requested this ``request`` fixture.
    '''

    # Parameter parametrizing the fixture requesting this "request" fixture
    fixture_param = None

    # If the caller passed a default value *AND* this fixture is unparametrized,
    # default this parameter to this value.
    if (default_value is not _DEFAULT_VALUE_NONE and
        not is_fixture_parametrized(request)):
        fixture_param = default_value
    # Else...
    else:
        # Raise an exception unless this fixture is parametrized.
        die_unless_fixture_parametrized(request)

        # Parameter parametrizing this fixture.
        fixture_param = request.param

    # Else, that fixture was parametrized. If the caller requested that
    # parametrization be validated, do so.
    if type_expected is not None and not isinstance(
        fixture_param, type_expected):
        raise BetseTestFixtureException(
            'Fixture "{}" parametrization not of type "{}": {!r}'.format(
                get_fixture_name(request), type_expected, fixture_param))

    # Return this parameter.
    return fixture_param

# ....................{ GETTERS ~ request                  }....................
#FIXME: This function internally calls the now-deprecated
#request.getfuncargvalue() function, which has since been replaced by the
#request.getfixturevalue() function. Sadly, the latter is unavailable in the
#most recent stable version of "pytest" available under my Linux distribution
#(Gentoo for the win) and hence ignored for the moment.
#
#Revisit this in early 2017, please.

@type_check
def get_requested_fixture(
    request,
    fixture_name: str,
    type_expected: TestableTypes = None,
) -> object:
    '''
    Object returned by the fixture with the passed name transitively requested
    by the current test requesting the passed ``request`` fixture object,
    optinionally validated to be an instance of the passed type or tuple of
    types.

    Parameters
    ----------
    request : _pytest.python.FixtureRequest
        Object passed to fixtures and tests requesting the ``request`` fixture.
    fixture_name: str
        Name of the fixture to return the object returned by that fixture.
    type_expected : optional[TestableTypes]
        Class or tuple of classes to validate this object to be an instance of.
        If this object is _not_ an instance of this class or classes, an
        exception is raised. Defaults to `None`, in which case this object is
        returned unvalidated.

    Returns
    ----------
    _pytest.python.FixtureDef
        Object returned by this fixture.

    Raises
    ----------
    `pytests.get_fixture_lookup_error_type()`
        If this fixture is either unavailable _or_ is available but
        unretrievable (e.g., this fixture is the child of the current fixture).
    '''

    # Object returned by this fixture. Terrible Function Names Part XIXI: the
    # endless journey endures.
    fixture_object = request.getfuncargvalue(fixture_name)

    # If the caller requested this object be validated, do so.
    if type_expected is not None and not isinstance(
        fixture_object, type_expected):
        raise BetseTestFixtureException(
            'Fixture "{}" object not of type "{}": {!r}'.format(
                fixture_name, type_expected, fixture_object))

    # Return this object.
    return fixture_object


def get_requested_fixture_or_none(
    request: '_pytest.python.FixtureRequest',
    fixture_name: str,
) -> '_pytest.python.FixtureDef':
    '''
    Fixture with the passed name transitively requested by the current test if
    this fixture is both available and retrievable _or_ `None` otherwise,
    inspected from the passed ``request`` fixture object.

    Parameters
    ----------
    request : _pytest.python.FixtureRequest
        Object passed to fixtures and tests requesting the ``request`` fixture.
    fixture_name: str
        Name of the fixture to return.

    Returns
    ----------
    _pytest.python.FixtureDef or None
        This fixture if this fixture is both available and retrievable _or_
        `None` otherwise.
    '''

    # Class of all "FixtureLookupError" exceptions raised by "pytest".
    FixtureLookupError = pytests.get_fixture_lookup_error_type()

    # Attempt to retrieve and return this fixture for this request.
    try:
        return get_requested_fixture(request, fixture_name)
    # If no such fixture exists, return None.
    except FixtureLookupError:
        return None

# ....................{ GETTERS ~ request : name           }....................
def get_requested_fixture_names(
    request: '_pytest.python.FixtureRequest') -> list:
    '''
    List of the names of all fixtures transitively requested by the current test
    excluding the fixture that explicitly requested the passed ``request`` fixture
    object if any, inspected from that ``request`` fixture object.

    This list includes the names of all fixtures that are:

    * Directly requested by this test.
    * Directly or indirectly requested by fixtures requested by this test.

    Parameters
    ----------
    request : _pytest.python.FixtureRequest
        Object passed to fixtures and tests requesting the ``request`` fixture.

    Returns
    ----------
    list
        List of all such names.
    '''

    # List of the names of all fixtures transitively required by this test,
    # including the fixture that explicitly requested this "request" fixture if
    # any. Since the original list is guaranteed *NOT* to be modified below,
    # this list need not be copied from the original list to preserve the
    # latter.
    #
    # This attribute is undocumented. Why is this attribute undocumented? Work
    # with me here, people.
    fixture_names = request.fixturenames

    # If a fixture explicitly requested this "request" fixture...
    if is_fixture(request):
        # Name of this fixture.
        omit_fixture_name = get_fixture_name(request)

        # Exclude this name from this list *WITHOUT* modifying the original
        # list, as doing so would subsequently raise "KeyError" exceptions.
        fixture_names = sequences.omit_item(fixture_names, omit_fixture_name)
        assert omit_fixture_name not in fixture_names

    # Return these fixture names.
    return fixture_names

# ....................{ GETTERS ~ request : name : prefix  }....................
def get_requested_fixture_name_prefixed_by(
    request: '_pytest.python.FixtureRequest',
    fixture_name_prefix: str,
) -> str:
    '''
    Name prefixed by the passed prefix of the single fixture transitively
    requested by the current test, inspected from the passed ``request`` fixture
    object.

    If either no such fixture _or_ more than one such fixture exist, an
    exception is raised.

    Parameters
    ----------
    request : _pytest.python.FixtureRequest
        Object passed to fixtures and tests requesting the ``request`` fixture.
    fixture_name_prefix: str
        String prefixing the fixture name to be returned.

    Returns
    ----------
    str
        Such name.

    See Also
    ----------
    :func:`get_requested_fixture_names`
        Further details.
    '''

    # List of the names of all fixtures prefixed by this prefix.
    prefixed_fixture_names = get_requested_fixture_names_prefixed_by(
        request, fixture_name_prefix)

    # Number of such fixtures.
    prefixed_fixture_count = len(prefixed_fixture_names)

    # If either no or more than one such fixtures exist, raise an exception.
    if prefixed_fixture_count != 1:
        # Exception message to be raised.
        exception_message = None

        if prefixed_fixture_count == 0:
            exception_message = (
                'No fixture prefixed by "{}" requested by this test.'.format(
                    fixture_name_prefix))
        else:
            exception_message = (
                'Multiple fixtures prefixed by "{}" '
                'requested by this test: {}'.format(
                    fixture_name_prefix, prefixed_fixture_names))

        # Raise this exception with this message.
        raise BetseTestFixtureException(exception_message)

    # Else, return the single such fixture name.
    return prefixed_fixture_names[0]


def get_requested_fixture_names_prefixed_by(
    request: '_pytest.python.FixtureRequest',
    fixture_name_prefix: str,
) -> list:
    '''
    List of all names prefixed by the passed prefix of all fixtures transitively
    requested by the current test, inspected from the passed ``request`` fixture
    object.

    Parameters
    ----------
    request : _pytest.python.FixtureRequest
        Object passed to fixtures and tests requesting the ``request`` fixture.
    fixture_name_prefix: str
        String prefixing all fixture names in the returned list.

    Returns
    ----------
    list
        List of all such names.

    See Also
    ----------
    get_requested_fixture_names
        Further details.
    '''

    return sequences.get_items_prefixed_by(
        sequence=get_requested_fixture_names(request),
        item_prefix=fixture_name_prefix,
    )
