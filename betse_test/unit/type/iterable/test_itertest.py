#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Unit tests exercising the :mod:`betse.util.type.iterable.itertest` submodule.
'''

# ....................{ IMPORTS                           }....................

# ....................{ TESTS ~ default                   }....................
#FIXME: Exercise the additional edge-case of a non-empty non-set of an
#arbitrary number of items, at least one of which (ideally the last) is
#non-hashable. Which begs the question, of course: what's trivially
#non-hashable?
def test_is_items_unique() -> None:
    '''
    Unit test the
    :func:`betse.util.type.iterable.itertest.is_items_unique` function.
    '''

    # Defer heavyweight imports.
    from betse.util.type.iterable import generators, itertest

    # Short iterable containing only unique items.
    mysidian_legend_unique = (
        'One to be born',
        'From a dragon',
        'Hoisting the light',
        'And the dark',
        'Arises high up',
        'In the sky to',
        'The still land.',
        'Veiling the moon with',
        'The Light of Eternity,',
        'It brings',
        'Another promise',
        'To Mother Earth with',
        'A bounty and mercy..',
    )

    # Short iterable containing one or more duplicate items.
    mysidian_legend_nonunique = (
        mysidian_legend_unique + (mysidian_legend_unique[-1],))

    # Long iterable containing only unique items.
    burning_legend_unique = tuple(range(256))

    # Long iterable containing one or more duplicate items.
    burning_legend_nonunique = (
        burning_legend_unique + (burning_legend_unique[-1],))

    # Exercise the edge case when this iterable is an empty generator.
    assert itertest.is_items_unique(iterable=generators.empty_generator())

    # Exercise the edge case when this iterable is an empty non-generator.
    assert itertest.is_items_unique(iterable=())

    # Exercise the edge case when this iterable is a non-empty set.
    assert itertest.is_items_unique(iterable=set(mysidian_legend_unique))

    # Exercise the edge case when this iterable is a non-empty non-set of less
    # than 64 hashable unique items.
    assert itertest.is_items_unique(iterable=mysidian_legend_unique)

    # Exercise the edge case when this iterable is a non-empty non-set of less
    # than 64 hashable non-unique items.
    assert (
        itertest.is_items_unique(iterable=mysidian_legend_nonunique) is False)

    # Exercise the edge case when this iterable is a non-empty non-set of
    # greater than or equal to than 64 hashable unique items.
    assert itertest.is_items_unique(iterable=burning_legend_unique)

    # Exercise the edge case when this iterable is a non-empty non-set of
    # greater than or equal to than 64 hashable non-unique items.
    assert (
        itertest.is_items_unique(iterable=burning_legend_nonunique) is False)
