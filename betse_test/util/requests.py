#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Utility functions simplifying access to py.test's built-in `request` fixture.

The official documentation for this fixture is blatantly inadequate. Even were
this documentation adequate, however, this fixture's API is overly obtuse to an
absurd (and almost obfuscatory) degree. This module provides an intelligible
alternative to that API.

See Also
----------
https://pytest.org/latest/builtin.html#_pytest.python.FixtureRequest
    Official documentation for this fixture – such as it is.
'''

# ....................{ IMPORTS                            }....................
import copy
from _pytest.python import FixtureLookupError
from betse.util.type import sequences, types
from betse_test.util.exceptions import BetseTestFixtureException

# ....................{ TESTERS                            }....................
def is_fixture(request: '_pytest.python.FixtureRequest') -> bool:
    '''
    `True` only if the the passed `request` fixture object was requested by a
    fixture rather than a test.

    Parameters
    ----------
    request : _pytest.python.FixtureRequest
        Object passed to fixtures and tests requesting the `request` fixture.

    Returns
    ----------
    bool
        `True` only if a fixture requested this `request` fixture.
    '''

    return request.fixturename is not None

# ....................{ GETTERS                            }....................
# All functions accepting a "request" fixture object assume the passed object to
# be of the expected type rather than validating this to be the case. Why?
# Because this type resides in a private module that does *NOT* appear to be
# externally importable: e.g.,
#
#     >>> import _pytest.python.FixtureRequest as req
#     ImportError: cannot import name 'transfer_markers'

def get_test_name(request: '_pytest.python.FixtureRequest') -> str:
    '''
    Unqualified name of the current test (e.g., `test_cli_info`), inspected from
    the passed `request` fixture object.

    This function returns either:

    * If the current test is a test function, this function's unqualified name.
    * If the current test is a test method, this method's unqualified name.

    Parameters
    ----------
    request : _pytest.python.FixtureRequest
        Object passed to fixtures and tests requesting the `request` fixture.

    Returns
    ----------
    str
        Unqualified name of this test.
    '''

    # No one would ever think of trying this. Except one man did.
    return request.node.name

# ....................{ GETTERS ~ fixture                  }....................
def get_fixture(
    request: '_pytest.python.FixtureRequest',
    fixture_name: str,
) -> '_pytest.python.FixtureDef':
    '''
    Fixture with the passed name transitively requested by the current test,
    inspected from the passed `request` fixture object.

    Parameters
    ----------
    request : _pytest.python.FixtureRequest
        Object passed to fixtures and tests requesting the `request` fixture.
    fixture_name: str
        Name of the fixture to return.

    Returns
    ----------
    _pytest.python.FixtureDef
        This fixture.

    Raises
    ----------
    _pytest.python.FixtureLookupError
        If this fixture is either unavailable _or_ is available but
        unretrievable (e.g., this fixture is the child of the current fixture).
    '''
    assert types.is_str_nonempty(fixture_name), (
        types.assert_not_str_nonempty(fixture_name, 'Fixture name'))

    # Terrible Function Names Part XIXI: the endless journey endures.
    return request.getfuncargvalue(fixture_name)


def get_fixture_or_none(
    request: '_pytest.python.FixtureRequest',
    fixture_name: str,
) -> '_pytest.python.FixtureDef':
    '''
    Fixture with the passed name transitively requested by the current test if
    this fixture is both available and retrievable _or_ `None` otherwise,
    inspected from the passed `request` fixture object.

    Parameters
    ----------
    request : _pytest.python.FixtureRequest
        Object passed to fixtures and tests requesting the `request` fixture.
    fixture_name: str
        Name of the fixture to return.

    Returns
    ----------
    _pytest.python.FixtureDef, None
        This fixture if this fixture is both available and retrievable _or_
        `None` otherwise..
    '''

    try:
        return get_fixture(request, fixture_name)
    except FixtureLookupError:
        return None

# ....................{ GETTERS ~ fixture : name           }....................
def get_fixture_name(request: '_pytest.python.FixtureRequest') -> str:
    '''
    Name of the current fixture that explicitly requested the passed `request`
    fixture object.

    If this `request` fixture was requested by a test rather than a fixture,
    an exception is raised.

    Parameters
    ----------
    request : _pytest.python.FixtureRequest
        Object passed to fixtures and tests requesting the `request` fixture.

    Returns
    ----------
    str
        Name of such fixture.
    '''

    # If this "request" fixture was requested by a test rather than fixture,
    # raise an exception.
    if not is_fixture(request):
        raise BetseTestFixtureException(
            '"request" fixture requested by '
            'a test rather than fixture: {}'.format(request))

    # Else, this "request" fixture was requested by a fixture. Return this
    # fixture's name.
    return request.fixturename


def get_fixture_name_prefixed_by(
    request: '_pytest.python.FixtureRequest',
    fixture_name_prefix: str,
) -> str:
    '''
    Name prefixed by the passed prefix of the single fixture transitively
    requested by the current test, inspected from the passed `request` fixture
    object.

    If either no such fixture or more than one such fixture exist, an exception
    is raised.

    Parameters
    ----------
    request : _pytest.python.FixtureRequest
        Object passed to fixtures and tests requesting the `request` fixture.
    fixture_name_prefix: str
        String prefixing the fixture name to be returned.

    Returns
    ----------
    str
        Such name.

    See Also
    ----------
    get_fixture_names
        Further details.
    '''

    # List of the names of all fixtures prefixed by this prefix.
    prefixed_fixture_names = get_fixture_names_prefixed_by(
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

# ....................{ GETTERS ~ fixture : names          }....................
def get_fixture_names(request: '_pytest.python.FixtureRequest') -> list:
    '''
    List of the names of all fixtures transitively requested by the current test
    excluding the fixture that explicitly requested the passed `request` fixture
    object if any, inspected from that `request` fixture object.

    This list includes the names of all fixtures that are:

    * Directly requested by this test.
    * Directly or indirectly requested by fixtures requested by this test.

    Parameters
    ----------
    request : _pytest.python.FixtureRequest
        Object passed to fixtures and tests requesting the `request` fixture.

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


def get_fixture_names_prefixed_by(
    request: '_pytest.python.FixtureRequest',
    fixture_name_prefix: str,
) -> list:
    '''
    List of all names prefixed by the passed prefix of all fixtures transitively
    requested by the current test, inspected from the passed `request` fixture
    object.

    Parameters
    ----------
    request : _pytest.python.FixtureRequest
        Object passed to fixtures and tests requesting the `request` fixture.
    fixture_name_prefix: str
        String prefixing all fixture names in the returned list.

    Returns
    ----------
    list
        List of all such names.

    See Also
    ----------
    get_fixture_names
        Further details.
    '''

    return sequences.get_items_prefixed_by(
        sequence=get_fixture_names(request),
        item_prefix=fixture_name_prefix,
    )
