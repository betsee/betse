#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Unit tests exercising the :mod:`betse.util.type.iterable.iterget` submodule.
'''

# ....................{ IMPORTS                           }....................

# ....................{ CLASSES                           }....................
class BellType(object):
    '''
    Type of each item of the iterable to be tested by the
    :func:`test_get_item_var_uniquified_str` unit test.

    Attributes
    ----------
    name : str
        Name of this item, guaranteed to be unique across all such items.
    stanza : str
        Arbitrary string associated with this item.
    '''

    # ..................{ INITIALIZERS                      }..................
    def __init__(self, name: str, stanza: str) -> None:

        self.name = name
        self.stanza = stanza

# ....................{ TESTS ~ default                   }....................
def test_get_item_var_uniquified_str() -> None:
    '''
    Unit test the
    :func:`betse.util.type.iterable.iterget.get_item_str_uniquified`
    function.
    '''

    # Defer heavyweight imports.
    from betse.util.type.iterable import iterget

    # Iterable to be tested.
    the_bells = [
        BellType(
            name='Silver bells !',
            stanza='What a world of merriment their melody foretells !'),
        BellType(
            name='Golden bells!',
            stanza='What a world of happiness their harmony foretells !'),
        BellType(
            name='Brazen bells !',
            stanza='What tale of terror, now, their turbulency tells !'),
        BellType(
            name='Iron bells !',
            stanza='What a world of solemn thought their monody compels !'),
    ]

    # Exercise the edge case when this iterable contains no item whose "name"
    # variable matches the passed format specifier.
    the_fifth_bell = iterget.get_item_str_uniquified(
        iterable=the_bells,
        item_attr_name='name',
        item_str_format='Runic bells ! ({})',
    )

    # Assert this getter function to have synthesized the expected string.
    assert the_fifth_bell == 'Runic bells ! (5)'

    # Add a new item whose "name" variable is this string to this iterable.
    the_bells.append(BellType(
        name=the_fifth_bell,
        stanza='In a clamorous appealing to the mercy of the fire,'))

    # Exercise the edge case when this iterable contains one item whose "name"
    # variable matches the passed format specifier.
    the_sixth_bell = iterget.get_item_str_uniquified(
        iterable=the_bells,
        item_attr_name='name',
        item_str_format='Runic bells ! ({})',
    )

    # Assert this getter function to have synthesized the expected string.
    assert the_sixth_bell == 'Runic bells ! (6)'
