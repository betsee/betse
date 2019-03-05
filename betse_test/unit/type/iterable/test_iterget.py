#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Unit tests exercising the :mod:`betse.util.type.iterable.iterget` submodule.
'''

# ....................{ IMPORTS                           }....................

# ....................{ CLASSES                           }....................
class RunicBells(object):
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
    :func:`betse.util.type.iterable.iterget.get_item_var_uniquified_str`
    function.
    '''

    # Defer heavyweight imports.
    from betse.util.type.iterable import iterget

    # Iterable to be tested.
    the_bells = (
        RunicBells(
            name='Silver bells !',
            stanza='What a world of merriment their melody foretells !'),
        RunicBells(
            name='Golden bells!',
            stanza='What a world of happiness their harmony foretells !'),
        RunicBells(
            name='Brazen bells !',
            stanza='What tale of terror, now, their turbulency tells !'),
        RunicBells(
            name='Iron bells !',
            stanza='What a world of solemn thought their monody compels !'),
    )

    # Exercise this getter function.
    the_fifth_bell = iterget.get_item_var_uniquified_str(
        iterable=the_bells,
        item_var_name='name',
        item_var_format='Runic bells ! ({})',
    )

    # Assert this getter function to have synthesized the expected string.
    assert the_fifth_bell == 'Runic bells ! (4)'
