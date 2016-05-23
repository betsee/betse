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
# from abc import ABCMeta

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
        Object passed to fixtures and tests requiring py.test's built-in
        `request` fixture.

    Returns
    ----------
    str
        Unqualified name of this test.
    '''

    # No one would ever think of trying this. Except one man did.
    return request.node.name


def get_test_fixture_names(request: '_pytest.python.FixtureRequest') -> list:
    '''
    List of the names of all fixtures transitively required by the current test,
    inspected from the passed `request` fixture object.

    This list includes the names of all fixtures that are:

    * Directly required by this test.
    * Directly or indirectly required by fixtures required by this test.

    Parameters
    ----------
    request : _pytest.python.FixtureRequest
        Object passed to fixtures and tests requiring py.test's built-in
        `request` fixture.

    Returns
    ----------
    list
        List of the names of all fixtures transitively required by this test.
    '''

    # This is undocumented. Why is this undocumented? Work with me here, people.
    return request.fixturenames
