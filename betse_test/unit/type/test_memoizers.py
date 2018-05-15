#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Unit tests exercising memoization decorators.
'''

# ....................{ IMPORTS                            }....................
# import pytest

# ....................{ TESTS                              }....................
def test_callable_cached() -> None:
    '''
    Unit test the :func:`betse.util.type.callables.func_cached` decorator.
    '''

    # Defer heavyweight imports.
    from betse.util.type.decorator.decmemo import func_cached

    # Class containing the callable to be cached.
    class Dreamland(object):
        def __init__(self, dreams: int) -> None:
            self.dreams = dreams

        @func_cached
        def dreams_cached(self) -> int:
            # Instance variable to be both cached and returned.
            dreams_cached = self.dreams * 2

            # To detect erroneous attempts to call this callable multiple
            # times for the same object, modify the instance variable producing
            # this return value in a predictable way on each call.
            self.dreams /= 2

            # Return this variable's value.
            return dreams_cached

    # Value of the "Keeper.dreams" attribute *BEFORE* invoking dreams_cached().
    DREAM_COUNT_PRECACHED = 7

    # Value of the "Keeper.dreams" attribute *AFTER* invoking dreams_cached().
    DREAM_COUNT_POSTCACHED = DREAM_COUNT_PRECACHED / 2

    # Value of the "Keeper.dreams_cached" callable.
    DREAM_COUNT_CACHED = DREAM_COUNT_PRECACHED * 2

    # Instance of this class initialized with this value.
    rebellion = Dreamland(dreams=DREAM_COUNT_PRECACHED)

    # Assert this attribute to be as initialized.
    assert rebellion.dreams == DREAM_COUNT_PRECACHED

    # Assert this callable to be as cached.
    assert rebellion.dreams_cached() == DREAM_COUNT_CACHED

    # Assert this attribute to have been modified by this callable call.
    assert rebellion.dreams == DREAM_COUNT_POSTCACHED

    # Assert this callable to still be as cached.
    assert rebellion.dreams_cached() == DREAM_COUNT_CACHED

    # Assert this attribute to *NOT* have been modified again, validating that
    # the prior callable access returned the previously cached value rather than
    # recalling this callable method.
    assert rebellion.dreams == DREAM_COUNT_POSTCACHED


def test_property_cached() -> None:
    '''
    Unit test the :func:`betse.util.type.callables.property_cached` decorator.
    '''

    # Defer heavyweight imports.
    from betse.util.type.decorator.decmemo import property_cached

    # Class containing the property to be cached.
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
