#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Unit tests exercising functionality defined by the
:mod:`betse.util.type.iterable.mapping.mapmerge` submodule.
'''

# ....................{ IMPORTS                           }....................
import pytest

# ....................{ TESTS                             }....................
def test_merge_maps() -> None:
    '''
    Unit test the :func:`betse.util.type.iterable.mapping.mapmerge.merge_maps`
    function.
    '''

    # Defer heavyweight imports.
    from betse.exceptions import BetseMappingException
    from betse.util.type.iterable.mapping import mapmerge
    from betse.util.type.iterable.mapping.mapmerge import MergeCollisionPolicy

    # Iterable of input dictionaries to be merged.
    hypnos = {
        'mercy':   'May the merciful gods, if indeed there be such,',
        'will':    'guard those hours when no power of the will,',
        'cunning': 'or drug that the cunning of man devises,',
        'sleep':   'can keep me from the chasm of sleep.',
    }
    the_call_of_cthulhu = {
        'mercy':     'The most merciful thing in the world, I think,',
        'inability': 'is the inability of the human mind',
        'correlate': 'to correlate all its contents.',
        'ignorance': 'We live on a placid island of ignorance',
        'infinity':  'in the midst of black seas of the infinity,',
        'voyage':    'and it was not meant that we should voyage far.',
    }
    the_colour_out_of_space = {
        'wild':  'West of Arkham the hills rise wild,',
        'woods': 'and there are valleys with deep woods',
        'cut':   'that no axe has ever cut.',

        # Key-value pair duplicated from the "hypnos" map to properly exercise
        # the "RAISE_EXCEPTION" merge policy.
        'mercy': 'May the merciful gods, if indeed there be such,',
    }
    the_other_gods = {
        'infinity': 'The vengeance of the infinite abysses...',
        'curse':    'That cursed, that damnable pit...',
        'mercy':    'Merciful gods of earth, I am falling into the sky!',
    }

    # Iterable of dictionaries containing one or more key collisions.
    mappings_nonunique = (hypnos, the_call_of_cthulhu, the_other_gods,)

    # Iterable of dictionaries containing *NO* key collisions.
    mappings_unique = (hypnos, the_colour_out_of_space)

    # Test whether merging dictionaries containing one or more key collisions
    # with exceptions raised on collisions raises the expected exception.
    with pytest.raises(BetseMappingException):
        mapmerge.merge_maps(
            mappings=mappings_nonunique,
            collision_policy=MergeCollisionPolicy.RAISE_EXCEPTION,
            is_values_copied=False,
        )

    # Test whether merging dictionaries containing no key collisions with
    # exceptions raised on collisions merges as expected.
    assert mapmerge.merge_maps(
        mappings=mappings_unique,
        collision_policy=MergeCollisionPolicy.RAISE_EXCEPTION,
        is_values_copied=True,
    ) == {
        'mercy':   'May the merciful gods, if indeed there be such,',
        'will':    'guard those hours when no power of the will,',
        'cunning': 'or drug that the cunning of man devises,',
        'sleep':   'can keep me from the chasm of sleep.',
        'wild':    'West of Arkham the hills rise wild,',
        'woods':   'and there are valleys with deep woods',
        'cut':     'that no axe has ever cut.',
    }

    # Test whether merging dictionaries containing key collisions with
    # precedence given to prior dictionaries succeeds as expected.
    assert mapmerge.merge_maps(
        mappings=mappings_nonunique,
        collision_policy=MergeCollisionPolicy.PREFER_FIRST,
        is_values_copied=False,
    ) == {
        'mercy':     'May the merciful gods, if indeed there be such,',
        'will':      'guard those hours when no power of the will,',
        'cunning':   'or drug that the cunning of man devises,',
        'sleep':     'can keep me from the chasm of sleep.',
        'inability': 'is the inability of the human mind',
        'correlate': 'to correlate all its contents.',
        'ignorance': 'We live on a placid island of ignorance',
        'infinity':  'in the midst of black seas of the infinity,',
        'voyage':    'and it was not meant that we should voyage far.',
        'curse':     'That cursed, that damnable pit...',
    }

    # Test whether merging dictionaries containing key collisions with
    # precedence given to subsequent dictionaries succeeds as expected.
    assert mapmerge.merge_maps(
        mappings=mappings_nonunique,
        collision_policy=MergeCollisionPolicy.PREFER_LAST,
        is_values_copied=True,
    ) == {
        'infinity':  'The vengeance of the infinite abysses...',
        'curse':     'That cursed, that damnable pit...',
        'mercy':     'Merciful gods of earth, I am falling into the sky!',
        'inability': 'is the inability of the human mind',
        'correlate': 'to correlate all its contents.',
        'ignorance': 'We live on a placid island of ignorance',
        'voyage':    'and it was not meant that we should voyage far.',
        'will':      'guard those hours when no power of the will,',
        'cunning':   'or drug that the cunning of man devises,',
        'sleep':     'can keep me from the chasm of sleep.',
    }
