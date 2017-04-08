#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Unit tests exercising the :class:`betse.util.type.mappings.OrderedArgsDict`
type.
'''

# ....................{ IMPORTS                            }....................
import pytest
from collections import OrderedDict

# ....................{ TESTS                              }....................
def test_orderedparamsdict_pass() -> None:
    '''
    Test all aspects of the `OrderedParmsDict` type intended to succeed.

    Parameters
    ----------
    betse_orderedparamsdict : OrderedArgsDict
        Instance of this type to be tested.
    '''

    # Defer heavyweight imports.
    from betse.util.type.mappings import OrderedArgsDict

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


def test_orderedparamsdict_fail() -> None:
    '''
    Test all aspects of the `OrderedArgsDict` type intended to fail.
    '''

    # Defer heavyweight imports.
    from betse.exceptions import BetseMappingException
    from betse.util.type.mappings import OrderedArgsDict

    # Test whether initializing this dictionary with an odd rather than even
    # number of positional arguments fails.
    with pytest.raises(BetseMappingException):
        OrderedArgsDict(
            'Danu', 'Don',
            'Dian Cecht',
        )
