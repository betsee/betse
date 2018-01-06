#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Unit tests exercising various simpler classes defined by the
:mod:`betse.util.type.mapping.mapcls` submodule *not* warranting their own
dedicated unit test submodule.
'''

# ....................{ IMPORTS                            }....................
from collections import OrderedDict
import pytest

# ....................{ TESTS ~ default                    }....................
def test_defaultdict_pass() -> None:
    '''
    Test aspects of the :class:`betse.util.type.mapping.mapcls.DefaultDict`
    type intended to succeed.
    '''

    # Defer heavyweight imports.
    from betse.util.type.mapping.mapcls import DefaultDict

    # Default dictionary to be tested.
    lake_isle_of_innisfree = DefaultDict(
        missing_key_value=lambda self, missing_key:
            missing_key + ' the deep heart’s core.',
        initial_mapping={
            'There midnight’s all': 'a glimmer',
            'and noon':             'a purple glow,',
            'And evening full':     'of the linnet’s wings.',
        },
    )

    # Test whether this dictionary contains a non-missing key initialized above.
    assert lake_isle_of_innisfree['and noon'] == 'a purple glow,'

    # Test whether a missing key defaults to the expected value.
    assert lake_isle_of_innisfree['I hear it in'] == (
        'I hear it in the deep heart’s core.')

# ....................{ TESTS ~ ordered                    }....................
def test_orderedargsdict_pass() -> None:
    '''
    Test aspects of the :class:`betse.util.type.mapping.mapcls.OrderedArgsDict`
    type intended to succeed.
    '''

    # Defer heavyweight imports.
    from betse.util.type.mapping.mapcls import OrderedArgsDict

    # Custom ordered dictionary to be tested.
    tuatha_de_danann = OrderedArgsDict(
        'Nuada', 'Nodens',
        'Lugh', 'Lugus',
        'Brigit', 'Brigantia',
    )

    # Standard ordered dictionary containing the same key-value pairs.
    tuath_de = OrderedDict((
        ('Nuada', 'Nodens'),
        ('Lugh', 'Lugus'),
        ('Brigit', 'Brigantia'),
    ))

    # Test whether these dictionaries contain the same key-value pairs.
    assert tuatha_de_danann == tuath_de

    # Test whether initializing an empty custom ordered dictionary succeeds.
    empty_dict = OrderedArgsDict()

    # Test whether this dictionary is indeed empty.
    assert not empty_dict


def test_orderedargsdict_fail() -> None:
    '''
    Test aspects of the :class:`betse.util.type.mapping.mapcls.OrderedArgsDict`
    type intended to fail.
    '''

    # Defer heavyweight imports.
    from betse.exceptions import BetseMappingException
    from betse.util.type.mapping.mapcls import OrderedArgsDict

    # Test whether initializing this dictionary with an odd rather than even
    # number of positional arguments fails.
    with pytest.raises(BetseMappingException):
        OrderedArgsDict(
            'Danu', 'Don',
            'Dian Cecht',
        )
