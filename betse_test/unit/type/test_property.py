#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Unit tests exercising object properties.
'''

# ....................{ IMPORTS                            }....................
# import pytest

# ....................{ TESTS                              }....................
def test_property_cached() -> None:
    '''
    Unit test the :func:`betse.util.type.objects.property_cached` decorator.
    '''

    # Defer heavyweight imports.
    from betse.util.type.obj.objs import property_cached

    # Unannotated function to be type checked.
    class Keeper(object):
        def __init__(self, keys: int) -> None:
            self.keys = keys

        @property_cached
        def keys_cached(self) -> int:
            # Property value to be both cached and returned.
            keys_cached = self.keys * 2

            # To detect erroneous attempts to call this property method multiple
            # times for the same object, modify the object variable from which
            # this property value derived in a predictable way on each call of
            # this property method.
            self.keys /= 2

            # Return this property value.
            return keys_cached

    # Value of the "Keeper.keys" attribute *BEFORE* invoking keys_cached().
    KEY_COUNT_PRECACHED = 7

    # Value of the "Keeper.keys" attribute *AFTER* invoking keys_cached().
    KEY_COUNT_POSTCACHED = KEY_COUNT_PRECACHED / 2

    # Value of the "Keeper.keys_cached" property.
    KEY_COUNT_CACHED = KEY_COUNT_PRECACHED * 2

    # Instance of this class initialized with this value.
    i_want_out = Keeper(keys=KEY_COUNT_PRECACHED)

    # Assert this attribute to be as initialized.
    assert i_want_out.keys == KEY_COUNT_PRECACHED

    # Assert this property to be as cached.
    assert i_want_out.keys_cached == KEY_COUNT_CACHED

    # Assert this attribute to have been modified by this property call.
    assert i_want_out.keys == KEY_COUNT_POSTCACHED

    # Assert this property to still be as cached.
    assert i_want_out.keys_cached == KEY_COUNT_CACHED

    # Assert this attribute to *NOT* have been modified again, validating that
    # the prior property access returned the previously cached value rather than
    # recalling this property method.
    assert i_want_out.keys == KEY_COUNT_POSTCACHED
