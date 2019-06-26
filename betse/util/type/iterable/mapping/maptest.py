#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level **mapping tester** (i.e., utility functions testing and validating
dictionary-like types or instances) functionality.
'''

# ....................{ IMPORTS                           }....................
from betse.exceptions import (
    BetseMappingException,
    BetseMappingKeyException,
    BetseMappingValueException,
)
from betse.util.type.types import (
    type_check, IterableTypes, MappingType, HashableType)

# ....................{ EXCEPTIONS ~ map : keys           }....................
@type_check
def die_unless_has_keys(mapping: MappingType, keys: IterableTypes) -> None:
    '''
    Raise an exception unless the passed dictionary contains *all* passed keys.

    Equivalently, this function raises an exception if this dictionary does
    *not* contain one or more passed keys.

    Parameters
    ----------
    mapping : MappingType
        Dictionary to be validated.
    keys : IterableTypes
        Iterable of all keys to be tested for.

    Raises
    ----------
    BetseMappingKeyException
        If this dictionary does *not* contain one or more passed keys.

    See Also
    ----------
    :func:`has_keys`
        Further details.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.text.string import strjoin

    # If this dictionary does *NOT* contain one or more passed keys...
    if not has_keys(mapping=mapping, keys=keys):
        # Set of all passed keys *NOT* in this dictionary.
        keys_missing = set(key for key in keys if key not in mapping)

        # Grammatically proper noun describing the number of such keys.
        keys_noun = 'key' if len(keys_missing) == 1 else 'keys'

        # Raise an exception embedding this set.
        raise BetseMappingKeyException(
            'Dictionary {} {} not found.'.format(
                keys_noun,
                strjoin.join_as_conjunction_double_quoted(*keys_missing)))

# ....................{ EXCEPTIONS ~ map : values         }....................
@type_check
def die_unless_values_unique(mapping: MappingType) -> None:
    '''
    Raise an exception unless all values of the passed dictionary are unique.

    Equivalently, this function raises an exception if any two key-value pairs
    of this dictionary share the same values.

    Parameters
    ----------
    mapping : MappingType
        Dictionary to be validated.

    Raises
    ----------
    BetseMappingValueException
        If at least one value of this dictionary is a duplicate.

    See Also
    ----------
    :func:`is_values_unique`
        Further details.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.iterable import iterget
    from betse.util.type.text.string import strjoin

    # If one or more values of this dictionary are duplicates...
    if not is_values_unique(mapping):
        # Set of all duplicate values in this dictionary.
        values_duplicate = iterget.get_items_duplicate(mapping.values())

        # Grammatically proper noun describing the number of such values.
        values_noun = 'value' if len(values_duplicate) == 1 else 'values'

        # Raise an exception embedding this set.
        raise BetseMappingValueException(
            'Dictionary {} {} duplicate.'.format(
                values_noun,
                strjoin.join_as_conjunction_double_quoted(*values_duplicate)))

# ....................{ EXCEPTIONS ~ maps                 }....................
@type_check
def die_if_maps_collide(*mappings: MappingType) -> None:
    '''
    Raise an exception if two or more passed dictionaries **collide** (i.e.,
    contain key-value pairs containing the same keys but differing values).

    Parameters
    ----------
    mappings : Tuple[MappingType]
        Tuple of all dictionaries to be validated.

    Raises
    ----------
    BetseMappingException
        If two or more passed dictionaries collide.

    See Also
    ----------
    :func:`is_maps_collide`
        Further details.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.iterable import iteriter
    from betse.util.type.text.string import strjoin

    # If two or more of these dictionaries collide...
    if is_maps_collide(*mappings):
        # Iterable of all possible pairs of these dictionaries.
        mappings_pairs = iteriter.iter_pairs(mappings)

        # For each possible pair of these dictionaries...
        for mappings_pair in mappings_pairs:
            # If this pair of dictionaries collides...
            if is_maps_collide(*mappings_pair):
                # Set of all key-value pairs unique to a single mapping.
                items_unique = mappings[0].items() ^ mappings[1].items()

                # Set of all keys visited while iterating this set.
                keys_visited = set()

                # Set of all non-unique keys shared by two or more such pairs.
                keys_collide = set()

                # For each key of such a pair...
                for key, _ in items_unique:
                    # If this key has already been visited, this is a
                    # non-unique key shared by two or more such pairs.
                    if key in keys_visited:
                        keys_collide.add(key)

                    # Mark this key as having been visited.
                    keys_visited.add(key)

                # Grammatically proper noun describing the number of such keys.
                keys_noun = 'key' if len(keys_collide) == 1 else 'keys'

                # Raise an exception embedding this set.
                raise BetseMappingException(
                    'Dictionary {} {} not unique.'.format(
                        keys_noun, strjoin.join_as_conjunction_double_quoted(
                            *keys_collide)))

# ....................{ EXCEPTIONS ~ maps : keys          }....................
@type_check
def die_unless_maps_keys_equal(*mappings: MappingType) -> None:
    '''
    Raise an exception unless all of the passed dictionaries contain the exact
    same keys.

    Equivalently, this function raises an exception if any key of any passed
    dictionary is *not* a key of any other such dictionary.

    Parameters
    ----------
    mappings : Tuple[MappingType]
        Tuple of all dictionaries to be validated.

    Raises
    ----------
    BetseMappingKeyException
        If any key of any passed dictionary is *not* a key of any other such
        dictionary.

    See Also
    ----------
    :func:`is_keys_equal`
        Further details.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.text.string import strjoin

    # If one or more of these dictionaries contain differing keys...
    if not is_maps_keys_equal(*mappings):
        # First passed mapping. Since the is_keys_equal() function necessarily
        # returns true if either no mappings or only one mapping are passed,
        # this function returning false implies that two or more mappings are
        # passed. Ergo, this mapping is guaranteed to exist.
        mapping_first = mappings[0]

        # For each mapping excluding the first...
        for mapping in mappings[1:]:
            # If the keys of this mapping differ from those of the first...
            if not is_maps_keys_equal(mapping, mapping_first):
                # Set of all keys differing between these two mappings.
                keys_unequal = mapping.keys().symmetric_difference(
                    mapping_first.keys())

                # Grammatically proper noun describing the number of such keys.
                keys_noun = 'key' if len(keys_unequal) == 1 else 'keys'

                # Raise an exception embedding this set.
                raise BetseMappingKeyException(
                    'Dictionary {} {} differ.'.format(
                        keys_noun, strjoin.join_as_conjunction_double_quoted(
                            *keys_unequal)))

# ....................{ TESTERS ~ map : keys              }....................
@type_check
def has_keys(mapping: MappingType, keys: IterableTypes) -> bool:
    '''
    ``True`` only if the passed dictionary contains *all* passed keys.

    Parameters
    ----------
    mapping : MappingType
        Dictionary to be tested.
    keys : IterableTypes
        Iterable of all keys to be tested for.

    Returns
    ----------
    bool
        ``True`` only if this dictionary contains *all* passed keys.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.iterable import itertest

    # If any key is unhashable and hence *NOT* a valid key, raise an exception.
    itertest.die_unless_items_instance_of(iterable=keys, cls=HashableType)

    # Return true only if...
    return (
        # If only one key is passed, optimize this common edge case with the
        # standard idiom for testing key containment.
        keys[0] in mapping if len(keys) == 1 else
        # Else, two or more keys are passed. In this case, fallback to a
        # general-case strategy testing key containment in a single one-liner.
        # And yes: this is ridiculously awesome.
        set(keys).issubset(mapping)
    )

# ....................{ TESTERS ~ map : values            }....................
@type_check
def is_values_unique(mapping: MappingType) -> bool:
    '''
    ``True`` only if all values of the passed dictionary are **unique** (i.e.,
    if *no* two key-value pairs of this dictionary share the same values).

    Parameters
    ----------
    mapping : MappingType
        Dictionary to be inspected.

    Returns
    ----------
    bool
        ``True`` only if *all* values of this dictionary are unique.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.iterable import itertest

    # Defer to an existing low-level tester for sanity.
    return itertest.is_items_unique(mapping.values())

# ....................{ TESTERS ~ maps                    }....................
@type_check
def is_maps_collide(*mappings: MappingType) -> bool:
    '''
    ``True`` only if two or more passed dictionaries **collide** (i.e., contain
    key-value pairs containing the same keys but differing values).

    Equivalently, this function returns ``True`` only if no key of a key-value
    pair of one mapping is also a key of a key-value pair of another mapping
    such that the two values of these two pairs differ.

    Parameters
    ----------
    mappings : Tuple[MappingType]
        Tuple of all dictionaries to be tested.

    Returns
    ----------
    bool
        ``True`` only if if two or more passed dictionaries collide.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.iterable.set import sets

    # If two mappings are passed, prematurely optimize this common case by
    # directly testing these two mappings.
    #
    # Dismantled, this is:
    #
    # * "mappings[0].items() | mappings[1].items()", the set of 2-tuples of all
    #   distinct key-value pairs in these mappings.
    # * "len(mappings[0].items() | mappings[1].items())", the number of all
    #   distinct key-value pairs in these mappings.
    # * "mappings[0].keys()  | mappings[1].keys()", the set of 2-tuples of all
    #   distinct keys in these mappings.
    # * "len(mappings[0].keys()  | mappings[1].keys())", the number of all
    #   distinct keys in these mappings.
    # * "len(...) != len(...)", true only if the number of all distinct
    #   key-value pairs in these mappings differs from the number of distinct
    #   keys in these mappings. If true, then the former is guaranteed to be
    #   strictly larger than the latter (i.e., more distinct key-value pairs
    #   than keys exist), in which case one or more distinct key-value pairs
    #   must necessarily share the same key and hence collide.
    #
    # For example:
    #
    #     # A pair of colliding mappings.
    #     >>> mappings = ({'a': 42, 'b': 24}, {'a': 35, 'c': 53})
    #     >>> items = mappings[0].items() | mappings[1].items()
    #     >>> keys = mappings[0].keys() | mappings[1].keys()
    #     >>> items
    #     {('a', 42), ('a', 35), ('b', 24), ('c', 53)}
    #     >>> keys
    #     {'a', 'b', 'c'}
    #     >>> len(items) > keys
    #     True
    #
    #     # A pair of non-colliding mappings.
    #     >>> mappings = ({'a': 42, 'b': 24}, {'a': 42, 'c': 53})
    #     >>> items = mappings[0].items() | mappings[1].items()
    #     >>> keys = mappings[0].keys() | mappings[1].keys()
    #     >>> items
    #     {('a', 42), ('b', 24), ('c', 53)}
    #     >>> keys
    #     {'a', 'b', 'c'}
    #     >>> len(items) > keys
    #     False
    if len(mappings) == 2:
        return (
            len(mappings[0].items() | mappings[1].items()) !=
            len(mappings[0].keys()  | mappings[1].keys())
        )
    # Else if either no mappings or only one mapping are passed, return true.
    # See the is_maps_keys_equal() function for discussion on this edge case.
    elif len(mappings) < 2:
        return True
    # Else, three or more mappings are passed. In this case, defer to a
    # general-purpose algorithm supporting arbitrarily many mappings.

    # Sets of distinct key-value pairs and keys in these mappings.
    mappings_items = sets.make_union(mapping.items() for mapping in mappings)
    mappings_keys  = sets.make_union(mapping.keys()  for mapping in mappings)

    # Return true only if the number of distinct key-value pairs in these
    # mapping differs from and hence is strictly greater than the number of
    # distinct keys in these mappings. See above for further details.
    return len(mappings_items) != len(mappings_keys)

# ....................{ TESTERS ~ maps                    }....................
@type_check
def is_maps_keys_equal(*mappings: MappingType) -> bool:
    '''
    ``True`` only if all passed dictionaries contain the exact same keys.

    Parameters
    ----------
    mappings : Tuple[MappingType]
        Tuple of all dictionaries to be tested.

    Returns
    ----------
    bool
        ``True`` only if these dictionaries contain the exact same keys.
    '''

    # If two mappings are passed, prematurely optimize this common case by
    # directly testing the keys contained in these two mappings for equality.
    if len(mappings) == 2:
        return mappings[0].keys() == mappings[1].keys()
    # Else if either no mappings or only one mapping are passed, return true.
    # Why? Because:
    #
    # * If no mappings are passed, this edge case is functionally equivalent to
    #   the edge case in which one or more empty mappings (i.e., mappings
    #   containing no key-value pairs) are passed. Since the key containers for
    #   empty mappings are themselves empty, these key containers contain the
    #   exact same keys -- namely, none.
    # * If only one mapping is passed, this mapping by definition contains the
    #   exact same keys as itself.
    elif len(mappings) < 2:
        return True
    # Else, three or more mappings are passed. In this case, defer to a general
    # case algorithm.

    # Keys of the first passed mapping.
    mapping_first_keys = mappings[0].keys()

    # Return true only if...
    return all(
        # This mapping contains the same keys as the first such mapping...
        mapping.keys() == mapping_first_keys
        # For each mapping excluding the first.
        for mapping in mappings[1:]
    )
