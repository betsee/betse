#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Unit tests exercising the :func:`betse.util.type.cls.descriptors.expr_alias`
descriptor.
'''

# ....................{ IMPORTS                            }....................
import pytest
from pytest import fixture

# ....................{ FIXTURES                           }....................
@fixture(scope='session')
def betse_expralias() -> object:
    '''
    Fixture creating and returning a mock object whose class contains instances
    of the :func:`betse.util.type.cls.descriptors.expr_alias` descriptor to be
    tested.
    '''

    # Defer heavyweight imports.
    from betse.util.type.cls.descriptors import expr_alias

    # Class containing instances of this descriptor.
    class SongOfAragorn(object):
        # Descriptor aliasing an existing dictionary entry with typing.
        line_four = expr_alias(
            expr='self._deep_roots["are not"]["reached by"]', cls=str)

        # Descriptor aliasing an existing dictionary entry with typing, whose
        # corresponding value is erroneously of a different type.
        line_two = expr_alias(
            expr='self._deep_roots["Not all those"]["who wander"]', cls=str)

        # Descriptor aliasing an existing dictionary entry without typing.
        line_five = expr_alias(
            expr='self._deep_roots["From the ashes"]["a fire shall"]')

        # Descriptor aliasing a non-existing dictionary entry with typing.
        line_none = expr_alias(
            expr='self._deep_roots["Many a man"]["his life hath"]', cls=str)

        # Initialize the nested dictionary structures aliased above.
        def __init__(self) -> None:
            self._deep_roots = {
                "are not": {
                    "reached by": "the frost.",
                },
                "Not all those": {
                    "who wander": 0xA3E01054,
                },
                "From the ashes": {
                    "a fire shall": "be woken,",
                }
            }

    # Create and return an instance of this type.
    return SongOfAragorn()

# ....................{ TESTS                              }....................
def test_expralias_pass(betse_expralias) -> None:
    '''
    Test all aspects of the :func:`betse.util.type.cls.descriptors.expr_alias`
    descriptor intended to succeed.

    Parameters
    ----------
    betse_expralias : object
        Instance of this type to be tested.
    '''

    # Test the typed data descriptor's __get__() implementation as an instance
    # rather than class variable.
    assert betse_expralias.line_four == "the frost."

    # Test the untyped data descriptor's __get__() implementation as an instance
    # rather than class variable.
    assert betse_expralias.line_five == "be woken,"

    # Test passing the typed data descriptor's __set__() implementation a valid
    # string type.
    line_two = "Not all those who wander are lost;"
    betse_expralias.line_four = line_two
    assert betse_expralias.line_four == line_two
    assert betse_expralias._deep_roots["are not"]["reached by"] == line_two

    # Test passing the untyped data descriptor's __set__() implementation an
    # arbitrary non-string type.
    betse_expralias.line_five = 0xCAFEBABE
    assert betse_expralias.line_five == 0xCAFEBABE
    assert betse_expralias._deep_roots["From the ashes"]["a fire shall"] == (
        0xCAFEBABE)


def test_expralias_fail(betse_expralias) -> None:
    '''
    Test all aspects of the :func:`betse.util.type.cls.descriptors.expr_alias`
    descriptor intended to fail.

    Parameters
    ----------
    betse_expralias : object
        Instance of this type to be tested.
    '''

    # Imports deferred for safety.
    from betse.exceptions import BetseTypeException

    # Test calling the typed data descriptor's __get__() implementation as a
    # class rather than instance variable. Since this implementation expects to
    # be called only as an instance variable, an exception is raised.
    with pytest.raises(AttributeError):
        betse_expralias.__class__.line_four

    # Test passing the typed data descriptor's __set__() implementation an
    # invalid non-string type.
    with pytest.raises(BetseTypeException):
        betse_expralias.line_four = 0xFEEDBABE

    # Test an invalid data descriptor's __get__() implementation, aliased to
    # access an existing key of an existing dictionary with an unexpected type.
    with pytest.raises(BetseTypeException):
        betse_expralias.line_two

    # Test an invalid data descriptor's __get__() implementation, aliased to
    # access a non-existing key of an existing dictionary.
    with pytest.raises(KeyError):
        betse_expralias.line_none

    # Test the invalid data descriptor's __set__() implementation, aliased to
    # assign a non-existing key of an existing dictionary.
    with pytest.raises(KeyError):
        betse_expralias.line_none = "Many a man his life hath sold"

    # Test calling an arbitrary data descriptor's __delete__() implementation,
    # which currently remains unimplemented to preserve backward compatibility.
    with pytest.raises(AttributeError):
        del betse_expralias.line_four
