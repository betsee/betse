#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level **mapping tester** (i.e., utility functions testing and validating
dictionary-like types or instances) functionality.
'''

# ....................{ IMPORTS                           }....................
from betse.exceptions import BetseMappingException
from betse.util.type.types import (
    type_check, IterableTypes, MappingType, HashableType)

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
        Dictionary to be inspected.

    Raises
    ----------
    BetseMappingException
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
        raise BetseMappingException(
            'Dictionary {} {} duplicate.'.format(
                values_noun,
                strjoin.join_as_conjunction_double_quoted(*values_duplicate)))

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
    BetseMappingException
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
                raise BetseMappingException(
                    'Dictionary {} {} differ.'.format(
                        keys_noun, strjoin.join_as_conjunction_double_quoted(
                            *keys_unequal)))


@type_check
def die_unless_maps_keys_unique(*mappings: MappingType) -> None:
    '''
    Raise an exception unless no passed dictionaries **collide** (i.e., contain
    the same key).

    Equivalently, this function raises an exception if any key of any passed
    dictionary is also a key of any other such dictionary.

    Parameters
    ----------
    mappings : Tuple[MappingType]
        Tuple of all dictionaries to be validated.

    Raises
    ----------
    BetseMappingException
        If any key of any passed dictionary is also a key of any other such
        dictionary.

    See Also
    ----------
    :func:`is_keys_unique`
        Further details.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.iterable import iteriter
    from betse.util.type.text.string import strjoin

    # If one or more of these dictionaries contain the same keys...
    if not is_maps_keys_unique(*mappings):
        # Iterable of all possible pairs of these dictionaries.
        mappings_pairs = iteriter.iter_pairs(mappings)

        # For each possible pair of these dictionaries...
        for mappings_pair in mappings_pairs:
            # If this pair of dictionaries contains the same keys...
            if not is_maps_keys_unique(*mappings_pair):
                # Set of all keys shared between these two mappings.
                keys_equal = mappings_pair[0].keys() & mappings_pair[1].keys()

                # Grammatically proper noun describing the number of such keys.
                keys_noun = 'key' if len(keys_equal) == 1 else 'keys'

                # Raise an exception embedding this set.
                raise BetseMappingException(
                    'Dictionary {} {} not unique.'.format(
                        keys_noun, strjoin.join_as_conjunction_double_quoted(
                            *keys_equal)))

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

# ....................{ TESTERS ~ maps : items            }....................

# ....................{ TESTERS ~ maps : keys             }....................
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


@type_check
def is_maps_keys_unique(*mappings: MappingType) -> bool:
    '''
    ``True`` only if no passed dictionaries **collide** (i.e., contain the same
    key).

    Equivalently, this function returns ``True`` only if all passed
    dictionaries contain unique keys.

    Parameters
    ----------
    mappings : Tuple[MappingType]
        Tuple of all dictionaries to be tested.

    Returns
    ----------
    bool
        ``True`` only if these dictionaries only contain unique keys.
    '''

    # If two mappings are passed, prematurely optimize this common case by
    # directly testing these two mappings for an empty key set intersection.
    if len(mappings) == 2:
        return not(mappings[0].keys() & mappings[1].keys())
    # Else if either no mappings or only one mapping are passed, return true.
    # See the is_keys_equal() function for discussion on this edge case.
    elif len(mappings) < 2:
        return True
    # Else, three or more mappings are passed. In this case, defer to a general
    # case algorithm.

    # Set of all keys of the first passed mapping. Note that the dict.keys()
    # object is *NOT* a set and hence does *NOT* provide the set.intersection()
    # method called below.
    mapping_keys_first = set(mappings[0].keys())

    # Generator comprehension iteratively yielding all keys of all remaining
    # passed mappings.
    mapping_keys_rest = (mapping.keys() for mapping in mappings[1:])

    # Return true only if the intersection of the set of all keys of the first
    # passed mapping with those of all subsequent mappings is the empty set.
    return not(mapping_keys_first.intersection(*mapping_keys_rest))
