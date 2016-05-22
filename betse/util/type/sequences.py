#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level **non-string sequence** (i.e., non-string object implementing the
abstract base class `collections.Sequence`) facilities.

See Also
----------
betse.util.type.types.is_sequence
    Further details on what constitutes sequences and non-string sequences.
'''

# ....................{ IMPORTS                            }....................
from betse.util.type import types

# ....................{ REPLACERS                          }....................
def replace_items(
    sequence: "collections.Sequence",
    item_replacements: "collections.Mapping") -> (
    'collections.Sequence'):
    '''
    Get a new non-string sequence transformed from the passed non-string
    sequence by replacing all elements of the latter equalling any key of the
    passed dictionary by that key's value.

    This method performs global search-and-replacement on non-string sequences.

    Alternatives
    ----------
    While list and tuple comprehensions provide a more efficient means of
    performing multiple global search-and-replacements on lists and tuples
    guaranteed to contain only strings, such comprehensions fail to generalize
    to arbitrary sequences and hence are avoided here. For posterity, this
    ancient technique resembles:

        >>> bestiary = ['lamia', 'succubus', 'sylph']
        >>> bestiarum = [
        ...     beast.replace('lamia', 'lamara').replace('sylph', 'nymph')
        ...     for beast in bestiary
        ... ]
        >>> print(bestiarum)
        ['lamara', 'succubus', 'nymph']

    Parameters
    ----------
    sequence : collections.Sequence
        Original sequence to be returned transformed. For generality, this
        sequence is _not_ modified.
    item_replacements : collections.Mapping
        Mapping whose:
        * Keys are input values to find in the passed sequence. For simplicity,
          keys are matched via object equality rather than more complex object
          matching alternatives (e.g., glob, regex, prefix, suffix, substring).
        * Values are output values to replace these input values with in the
          returned sequence.

    Returns
    ----------
    collections.Sequence
        Sequence transformed from the passed sequence. For efficiency, this
        sequence is only a shallow rather than deep copy of the passed sequence.
    '''
    assert types.is_sequence_nonstr(sequence), (
        types.assert_not_sequence_nonstr(sequence))
    assert types.is_mapping(item_replacements), (
        types.assert_not_mapping(item_replacements))

    # Type of both the passed sequence and the sequence to be returned.
    sequence_type = type(sequence)

    # Return a shallow copy of this sequence, replacing each element that is a
    # key of the passed mapping by that key's value.
    return sequence_type(
        item_replacements[item] if item in item_replacements else item
        for item in sequence
    )

# ....................{ SORTERS                            }....................
def sort_lexicographic_ascending(
    sequence: "collections.Sequence") -> 'collections.Sequence':
    '''
    Get a new non-string sequence sorted from the passed non-string sequence in
    **ascending lexicographic order** (i.e., traditional order of dead-tree
    dictionaries and encyclopedias).

    Parameters
    ----------
    sequence : collections.Sequence
        Unsorted sequence to be returned sorted. For generality, this sequence
        is _not_ modified.

    Returns
    ----------
    collections.Sequence
        Sequence sorted from the passed sequence. For efficiency, this sequence
        is only a shallow rather than deep copy of the passed sequence.
    '''
    assert types.is_sequence_nonstr(sequence), (
        types.assert_not_sequence_nonstr(sequence))

    # Well, that was easy.
    return sorted(sequence)
