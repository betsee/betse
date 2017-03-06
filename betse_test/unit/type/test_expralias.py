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
from enum import Enum
from pytest import fixture

# ....................{ GLOBALS                            }....................
# Enumeration type to be aliased by the subsequent class. As required by the
# expr_enum_alias() data descriptor, the names of all members of this
# enumeration *MUST* be uppercase.
_ElendilsOath = Enum('_ElendilsOath', (
    'SINOME', 'MARUVAN', 'AR', 'HILDINYAR', 'TENN' 'AMBAR', 'METTA',))

# Enumeration type failing to satisfy the requirements of the
# expr_enum_alias() data descriptor. Specifically, at least one of the names of
# a member of this enumeration is *NOT* uppercase.
_Namarie = Enum('_Namarie', (
    'NU', 'LUINI', 'YASSEN', 'TINTILAR', 'i', 'ELENI',))

# ....................{ FIXTURES                           }....................
@fixture(scope='session')
def betse_expralias() -> object:
    '''
    Fixture creating and returning a mock object whose class contains instances
    of all expression alias descriptors to be unit-tested.
    '''

    # Imports deferred for safety.
    from betse.util.type.cls.expralias import (
        ExprAliasBound, expr_alias, expr_enum_alias)

    # Class containing instances of this descriptor.
    class SongOfAragorn(object):
        # Descriptor aliasing an existing dictionary entry with class typing.
        line_four = expr_alias(
            expr='self._deep_roots["are not"]["reached by"]', cls=str)

        # Descriptor aliasing an existing dictionary entry with predicate
        # typing.
        line_one = expr_alias(
            expr='self._deep_roots["All that is"]["gold does not"]',
            predicate=lambda value:
                value in ('glamour,', 'glimmer,','glitter,', 'glisten,',),
            predicate_label='adjectival',
        )

        # Descriptor aliasing an existing dictionary entry with both class and
        # predicate expression typing.
        line_six = expr_alias(
            expr='self._deep_roots["A light from"]["the shadows shall"]',
            cls=str,
            predicate_expr=(
                "value in ('spring;', 'summer;', 'autumn;', 'winter;',)"),
            predicate_label='seasonal',
        )

        # Descriptor aliasing an existing dictionary entry without typing.
        line_five = expr_alias(
            expr='self._deep_roots["From the ashes"]["a fire shall"]')

        # Descriptor aliasing an existing dictionary entry with enumerability,
        # whose corresponding value is defined below to be a lowercase string
        # whose uppercased equivalent is the name of an arbitrary member of the
        # enumeration defined above.
        line_three = expr_enum_alias(
            expr='self._deep_roots["The old that is strong"]["does not wither"]',
            enum_type=_ElendilsOath)

        # Descriptor aliasing an existing dictionary entry with typing, whose
        # corresponding value is defined below to be of a different type.
        line_two = expr_alias(
            expr='self._deep_roots["Not all those"]["who wander"]', cls=str)

        # Descriptor aliasing an existing dictionary entry with enumerability,
        # whose corresponding value is defined below to be of a non-enumeration
        # member.
        line_seven = expr_enum_alias(
            expr='self._deep_roots["Renewed shall be"]["blade that was"]',
            enum_type=_ElendilsOath)

        # Descriptor aliasing a non-existing dictionary entry with typing.
        line_none = expr_alias(
            expr='self._deep_roots["Many a man"]["his life hath"]', cls=str)


        def __init__(self) -> None:

            # Initialize the nested dictionary structures aliased above.
            self._deep_roots = {
                "are not": {
                    "reached by": "the frost.",
                },
                "All that is": {
                    "gold does not": "glitter,",
                },
                "Not all those": {
                    "who wander": 0xA3E01054,
                },
                "The old that is strong": {
                    "does not wither": _ElendilsOath.HILDINYAR.name.lower(),
                },
                "From the ashes": {
                    "a fire shall": "be woken,",
                },
                "A light from": {
                    "the shadows shall": "spring;",
                },
                "Renewed shall be": {
                    "blade that was": "broken,",
                },
                "The crownless again": {
                    "shall be": "king.",
                },
            }

            # Alias an arbitrary entry of this dictionary for use below.
            the_crownless_again = self._deep_roots['The crownless again']

            # Object aliasing an existing dictionary entry with class typing.
            # For coverage, the optional "obj_name" parameter is also exercised.
            self.line_eight = ExprAliasBound(
                expr='THE_CROWNLESS_AGAIN["shall be"]',
                cls=str,
                obj=the_crownless_again,
                obj_name='THE_CROWNLESS_AGAIN',
            )

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

    # Test the class-typed data descriptor's __get__() method.
    assert betse_expralias.line_four == "the frost."

    # Test the predicate-typed data descriptor's __get__() method.
    assert betse_expralias.line_one == "glitter,"

    # Test the class- and predicate-typed data descriptor's __get__()
    # method.
    assert betse_expralias.line_six == "spring;"

    # Test the enumeration-typed data descriptor's __get__() method.
    assert betse_expralias.line_three is _ElendilsOath.HILDINYAR

    # Test the untyped data descriptor's __get__() method.
    assert betse_expralias.line_five == "be woken,"

    # Test the class-typed object's get() method.
    assert betse_expralias.line_eight.get() == "king."

    # Test passing the class-typed data descriptor's __set__() method another
    # string value.
    merchant_of_venice = "Gilded tombs do worms enfold."
    betse_expralias.line_four = merchant_of_venice
    assert betse_expralias.line_four == merchant_of_venice
    assert betse_expralias._deep_roots["are not"]["reached by"] == (
        merchant_of_venice)

    # Test passing the predicate-typed data descriptor's __set__() method
    # another string value.
    betse_expralias.line_one              = 'glisten,'
    assert betse_expralias.line_one      == 'glisten,'
    assert betse_expralias._deep_roots[
        "All that is"]["gold does not"] == ('glisten,')

    # Test passing the class- and predicate-typed data descriptor's __set__()
    # method another string value.
    betse_expralias.line_six                   = 'winter;'
    assert betse_expralias.line_six           == 'winter;'
    assert betse_expralias._deep_roots[
        "A light from"]["the shadows shall"] == ('winter;')

    # Test passing the enumeration-typed data descriptor's __set__() method a
    # valid enumeration member.
    betse_expralias.line_three = _ElendilsOath.MARUVAN
    assert betse_expralias.line_three == _ElendilsOath.MARUVAN
    assert betse_expralias._deep_roots[
        "The old that is strong"]["does not wither"] == (
        _ElendilsOath.MARUVAN.name.lower())

    # Test passing the untyped data descriptor's __set__() method an arbitrary
    # non-string type.
    betse_expralias.line_five = 0xCAFEBABE
    assert betse_expralias.line_five == 0xCAFEBABE
    assert betse_expralias._deep_roots["From the ashes"]["a fire shall"] == (
        0xCAFEBABE)

    # Test passing the class-typed object's set() method another string value.
    prince_of_morocco = "Many a man his life hath sold"
    betse_expralias.line_eight.set(prince_of_morocco)
    assert betse_expralias.line_eight.get() == prince_of_morocco
    assert betse_expralias._deep_roots["The crownless again"]["shall be"] == (
        prince_of_morocco)


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
    from betse.exceptions import BetseEnumException, BetseTypeException
    from betse.util.type.cls.expralias import expr_enum_alias

    # Test instantiating the enumeration-typed data descriptor with an
    # enumeration type failing to satisfy this descriptor's requirements.
    with pytest.raises(BetseEnumException):
        expr_enum_alias(
            expr='wherein the stars tremble in the song of her voice,',
            enum_type=_Namarie)

    # Test calling the class-typed data descriptor's __get__() implementation as
    # a class rather than instance variable. Since this implementation expects
    # to be called only as an instance variable, an exception is raised.
    with pytest.raises(AttributeError):
        betse_expralias.__class__.line_four

    # Test passing the class-typed data descriptor's __set__() implementation a
    # non-string value.
    with pytest.raises(BetseTypeException):
        betse_expralias.line_four = 0xFEEDBABE

    # Test passing the predicate-typed data descriptor's __set__()
    # implementation an unrecognized string value.
    with pytest.raises(BetseTypeException):
        betse_expralias.line_one = 'expurgate,'

    # Test passing the class- and -predicate-typed data descriptor's __set__()
    # implementation both a non-string value and an unrecognized string value.
    with pytest.raises(BetseTypeException):
        betse_expralias.line_six = 0xFEEDFACE
    with pytest.raises(BetseTypeException):
        betse_expralias.line_six = 'extenuate;'

    # Test passing the enumeration-typed data descriptor's __set__()
    # implementation a non-enumeration member.
    with pytest.raises(BetseEnumException):
        betse_expralias.line_three = _Namarie.TINTILAR

    # Test an invalid class-typed data descriptor's __get__() implementation,
    # aliased to access an existing key of an existing dictionary whose value
    # has an unexpected type.
    with pytest.raises(BetseTypeException):
        betse_expralias.line_two

    # Test an invalid enumeration-typed data descriptor's __get__()
    # implementation, aliased to access an existing key of an existing
    # dictionary whose value is *NOT* the lowercase name of an enumeration
    # member.
    with pytest.raises(BetseEnumException):
        betse_expralias.line_seven

    # Test an invalid untyped data descriptor's __get__() implementation,
    # aliased to access a non-existing key of an existing dictionary.
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
