#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Unit tests exercising the :func:`betse.util.type.cls.descriptors.expr_alias`
descriptor.
'''

# ....................{ IMPORTS                            }....................
import pytest
from enum import Enum
from pytest import fixture

# ....................{ SUPERCLASSES                       }....................
class ExprAliasBaseClass(object):
    '''
    Placeholder superclass of an arbitrary expression alias instantiated below.
    '''

    pass

# ....................{ GLOBALS                            }....................
_ElendilsOath = Enum('_ElendilsOath', (
    'SINOME', 'MARUVAN', 'AR', 'HILDINYAR', 'TENN' 'AMBAR', 'METTA',))
'''
Enumeration type to be aliased by the subsequent class. As required by the
:func:`expr_enum_alias` data descriptor, the names of all members of this
enumeration *must* be uppercase.
'''


_Namarie = Enum('_Namarie', (
    'NU', 'LUINI', 'YASSEN', 'TINTILAR', 'i', 'ELENI',))
'''
Enumeration type failing to satisfy the requirements of the
:func:`expr_enum_alias(`data descriptor. Specifically, at least one of the names
of a member of this enumeration is *not* uppercase.
'''

# ....................{ FIXTURES                           }....................
@fixture(scope='session')
def betse_expralias() -> object:
    '''
    Fixture creating and returning a mock object whose ad-hoc class contains
    instances of all expression alias descriptors to be unit tested.
    '''

    # Imports deferred for safety.
    from betse.util.type.descriptor.expralias import (
        expr_alias, expr_enum_alias)
    from betse.util.type.descriptor.datadescs import DataDescriptorBound

    # Class containing instances of this descriptor.
    class SongOfAragorn(object):
        # Descriptor aliasing an existing dictionary entry with class typing
        # and subclassing an ad-hoc base class.
        soa_line_four = expr_alias(
            expr='self._deep_roots["are not"]["reached by"]',
            cls=str,
            base_classes=(ExprAliasBaseClass,),
        )

        # Descriptor aliasing an existing dictionary entry with predicate
        # typing.
        soa_line_one = expr_alias(
            expr='self._deep_roots["All that is"]["gold does not"]',
            predicate=lambda value:
                value in ('glamour,', 'glimmer,','glitter,', 'glisten,',),
            predicate_label='adjectival',
        )

        # Descriptor aliasing an existing dictionary entry with both class and
        # predicate expression typing.
        soa_line_six = expr_alias(
            expr='self._deep_roots["A light from"]["the shadows shall"]',
            cls=str,
            predicate_expr=(
                "value in ('spring;', 'summer;', 'autumn;', 'winter;',)"),
            predicate_label='seasonal',
        )

        # Descriptor aliasing an existing dictionary entry without typing.
        soa_line_five = expr_alias(
            expr='self._deep_roots["From the ashes"]["a fire shall"]')

        # Descriptor aliasing an existing dictionary entry with enumerability,
        # whose corresponding value is defined below to be a lowercase string
        # whose uppercased equivalent is the name of an arbitrary member of the
        # enumeration defined above.
        soa_line_three = expr_enum_alias(
            expr='self._deep_roots["The old that is strong"]["does not wither"]',
            enum_type=_ElendilsOath)

        # Descriptor aliasing an existing dictionary entry with typing, whose
        # corresponding value is defined below to be of a different type safely
        # castable to the type declared here.
        or_line_three = expr_alias(
            expr='self._one_ring["Nine for Mortal Men"]["doomed to die,"]',
            cls=float)

        # Descriptor aliasing an existing dictionary entry with typing, whose
        # corresponding value is defined below to be of a different type *NOT*
        # safely castable to the type declared here.
        soa_line_two = expr_alias(
            expr='self._deep_roots["Not all those"]["who wander"]', cls=str)

        # Descriptor aliasing an existing dictionary entry with enumerability,
        # whose corresponding value is defined below to be of a non-enumeration
        # member.
        soa_line_seven = expr_enum_alias(
            expr='self._deep_roots["Renewed shall be"]["blade that was"]',
            enum_type=_ElendilsOath)

        # Descriptor aliasing a non-existing dictionary entry with typing.
        soa_line_none = expr_alias(
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

            self._one_ring = {
                "Nine for Mortal Men": {
                    "doomed to die,": 1,
                },
            }

            # Alias an arbitrary entry of this dictionary for use below.
            the_crownless_again = self._deep_roots['The crownless again']

            # Object aliasing an existing dictionary entry with class typing.
            # For coverage, the optional "obj_name" parameter is also exercised.
            self.soa_line_eight = DataDescriptorBound(
                obj=the_crownless_again,
                data_desc=expr_alias(
                    expr='THE_CROWNLESS_AGAIN["shall be"]',
                    cls=str,
                    obj_name='THE_CROWNLESS_AGAIN',
                ))

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
    assert betse_expralias.soa_line_four == "the frost."

    # Test the numeric class-typed data descriptor's __get__() method,
    # implicitly casting a value of a different type to the desired type.
    assert type(betse_expralias.or_line_three) is float
    assert betse_expralias.or_line_three == 1.0

    # Test the numeric class-typed data descriptor's "expr_alias_cls" class
    # variable, accessed from this class rather than this instance.
    assert betse_expralias.__class__.or_line_three.expr_alias_cls is float

    # Test the predicate-typed data descriptor's __get__() method.
    assert betse_expralias.soa_line_one == "glitter,"

    # Test the class- and predicate-typed data descriptor's __get__()
    # method.
    assert betse_expralias.soa_line_six == "spring;"

    # Test the enumeration-typed data descriptor's __get__() method.
    assert betse_expralias.soa_line_three is _ElendilsOath.HILDINYAR

    # Test the untyped data descriptor's __get__() method.
    assert betse_expralias.soa_line_five == "be woken,"

    # Test the class-typed object's get() method.
    assert betse_expralias.soa_line_eight.get() == "king."

    # Test passing the class-typed data descriptor's __set__() method another
    # value of the same type.
    merchant_of_venice = "Gilded tombs do worms enfold."
    betse_expralias.soa_line_four = merchant_of_venice
    assert betse_expralias.soa_line_four == merchant_of_venice
    assert betse_expralias._deep_roots["are not"]["reached by"] == (
        merchant_of_venice)

    # Test passing the numeric class-typed data descriptor's __set__() method a
    # value of a different type safely castable to the desired type.
    betse_expralias.or_line_three = 9
    assert type(betse_expralias.or_line_three) is float
    assert betse_expralias.or_line_three == 9.0
    assert betse_expralias._one_ring[
        "Nine for Mortal Men"]["doomed to die,"] == 9.0

    # Test passing the predicate-typed data descriptor's __set__() method
    # another value of the same type.
    betse_expralias.soa_line_one              = 'glisten,'
    assert betse_expralias.soa_line_one      == 'glisten,'
    assert betse_expralias._deep_roots[
        "All that is"]["gold does not"] == ('glisten,')

    # Test passing the class- and predicate-typed data descriptor's __set__()
    # method another value of the same type.
    betse_expralias.soa_line_six                   = 'winter;'
    assert betse_expralias.soa_line_six           == 'winter;'
    assert betse_expralias._deep_roots[
        "A light from"]["the shadows shall"] == ('winter;')

    # Test passing the enumeration-typed data descriptor's __set__() method a
    # valid enumeration member.
    betse_expralias.soa_line_three = _ElendilsOath.MARUVAN
    assert betse_expralias.soa_line_three == _ElendilsOath.MARUVAN
    assert betse_expralias._deep_roots[
        "The old that is strong"]["does not wither"] == (
        _ElendilsOath.MARUVAN.name.lower())

    # Test passing the untyped data descriptor's __set__() method an arbitrary
    # non-string type.
    betse_expralias.soa_line_five = 0xCAFEBABE
    assert betse_expralias.soa_line_five == 0xCAFEBABE
    assert betse_expralias._deep_roots["From the ashes"]["a fire shall"] == (
        0xCAFEBABE)

    # Test passing the class-typed object's set() method another string value.
    prince_of_morocco = "Many a man his life hath sold"
    betse_expralias.soa_line_eight.set(prince_of_morocco)
    assert betse_expralias.soa_line_eight.get() == prince_of_morocco
    assert betse_expralias._deep_roots["The crownless again"]["shall be"] == (
        prince_of_morocco)

    # Test calling the class-typed data descriptor's __get__() implementation as
    # a class rather than instance variable.
    assert isinstance(betse_expralias.__class__.soa_line_four, ExprAliasBaseClass)


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
    from betse.exceptions import (
        BetseExprAliasException, BetseEnumException, BetseTypeException)
    from betse.util.type.descriptor.expralias import expr_enum_alias

    # Test instantiating the enumeration-typed data descriptor with an
    # enumeration type failing to satisfy this descriptor's requirements.
    with pytest.raises(BetseEnumException):
        expr_enum_alias(
            expr='wherein the stars tremble in the song of her voice,',
            enum_type=_Namarie)

    # Test passing the class-typed data descriptor's __set__() implementation a
    # non-string value.
    with pytest.raises(BetseTypeException):
        betse_expralias.soa_line_four = 0xFEEDBABE

    # Test passing the predicate-typed data descriptor's __set__()
    # implementation an unrecognized string value.
    with pytest.raises(BetseTypeException):
        betse_expralias.soa_line_one = 'expurgate,'

    # Test passing the class- and -predicate-typed data descriptor's __set__()
    # implementation both a non-string value and an unrecognized string value.
    with pytest.raises(BetseTypeException):
        betse_expralias.soa_line_six = 0xFEEDFACE
    with pytest.raises(BetseTypeException):
        betse_expralias.soa_line_six = 'extenuate;'

    # Test passing the enumeration-typed data descriptor's __set__()
    # implementation a non-enumeration member.
    with pytest.raises(BetseEnumException):
        betse_expralias.soa_line_three = _Namarie.TINTILAR

    # Test an invalid class-typed data descriptor's __get__() implementation,
    # aliased to access an existing key of an existing dictionary whose value
    # has an unexpected type.
    with pytest.raises(BetseTypeException):
        betse_expralias.soa_line_two

    # Test an invalid enumeration-typed data descriptor's __get__()
    # implementation, aliased to access an existing key of an existing
    # dictionary whose value is *NOT* the lowercase name of an enumeration
    # member.
    with pytest.raises(BetseEnumException):
        betse_expralias.soa_line_seven

    # Test an invalid untyped data descriptor's __get__() implementation,
    # aliased to access a non-existing key of an existing dictionary.
    with pytest.raises(BetseExprAliasException):
        betse_expralias.soa_line_none

    # Test the invalid data descriptor's __set__() implementation, aliased to
    # assign a non-existing key of an existing dictionary.
    with pytest.raises(KeyError):
        betse_expralias.soa_line_none = "Many a man his life hath sold"

    # Test calling an arbitrary data descriptor's __delete__() implementation,
    # which currently remains unimplemented to preserve backward compatibility.
    with pytest.raises(AttributeError):
        del betse_expralias.soa_line_four
