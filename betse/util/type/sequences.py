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

# ....................{ GETTERS                            }....................
def get_items_satisfying(
    sequence: "collections.Sequence", item_satisfier: callable) -> (
    'collections.Sequence'):
    '''
    Get a new non-string sequence containing only the proper subset of elements
    from the passed non-string sequence satisfying the passed callable.

    This method effectively performs a general-purpose global search for
    sequence elements.

    Parameters
    ----------
    sequence : collections.Sequence
        Original sequence to return a proper subset of. For safety, this
        function does _not_ modify this sequence.
    item_satisfier : callable
        Callable (e.g., function, lambda) accepting a single element of this
        sequence and returning only:
        * `True` if this element satisfies the desired requirements.
        * `False` otherwise.

    Returns
    ----------
    collections.Sequence
        Proper subset of the passed sequence. For efficiency, this new sequence
        is only a shallow rather than deep copy of the passed sequence.
    '''
    assert types.is_sequence_nonstr(sequence), (
        types.assert_not_sequence_nonstr(sequence))
    assert types.is_callable(item_satisfier), (
        types.assert_not_callable(item_satisfier))

    # Type of both the passed sequence and the sequence to be returned.
    sequence_type = type(sequence)

    # Return a generator-based shallow copy of this sequence.
    return sequence_type(
        item
        for item in sequence
        if item_satisfier(item)
    )


# ....................{ GETTERS ~ str                      }....................
def get_items_prefixed_by(
    sequence: "collections.Sequence", item_prefix: str) -> (
    'collections.Sequence'):
    '''
    Get a new non-string sequence containing only the proper subset of elements
    from the passed non-string sequence that are strings prefixed by the passed
    string prefix.

    This method effectively performs a prefix-based global search for sequence
    elements.

    Parameters
    ----------
    sequence : collections.Sequence
        Original sequence to return a proper subset of. For safety, this
        function does _not_ modify this sequence.
    item_prefix : str
        String prefixing all elements of the returned sequence.

    Returns
    ----------
    collections.Sequence
        Proper subset of the passed sequence. For efficiency, this new sequence
        is only a shallow rather than deep copy of the passed sequence.
    '''
    assert types.is_str_nonempty(item_prefix), (
        types.assert_not_str_nonempty(item_prefix, 'Element prefix'))

    # Return a generator-based shallow copy of this sequence.
    return get_items_satisfying(
        sequence=sequence,
        item_satisfier=lambda item: (
            types.is_str(item) and item.startswith(item_prefix)
        ),
    )

# ....................{ OMITTERS                           }....................
def omit_item(
    sequence: "collections.Sequence", omit_item: object) -> (
    'collections.Sequence'):
    '''
    Get a new non-string sequence containing only the proper subset of elements
    from the first passed non-string sequence that do _not_ equal the passed
    object.

    Parameters
    ----------
    sequence : collections.Sequence
        Original sequence to return a proper subset of. For safety, this
        function does _not_ modify this sequence.
    omit_item : object
        Object to be omitted from the returned sequence.

    Returns
    ----------
    collections.Sequence
        Proper subset of the passed sequence. For efficiency, this new sequence
        is only a shallow rather than deep copy of the passed sequence.
    '''

    return omit_items(
        sequence=sequence,
        omit_items=tuple(omit_item),
    )


def omit_items(
    sequence: "collections.Sequence",
    omit_items: "collections.Sequence") -> (
    'collections.Sequence'):
    '''
    Get a new non-string sequence containing only the proper subset of elements
    from the first passed non-string sequence that do _not_ equal elements of
    the second passed non-string sequence.

    Parameters
    ----------
    sequence : collections.Sequence
        Original sequence to return a proper subset of. For safety, this
        function does _not_ modify this sequence.
    omit_items : collections.Sequences
        Sequence containing all elements to be omitted from the returned
        sequence.

    Returns
    ----------
    collections.Sequence
        Proper subset of the passed sequence. For efficiency, this new sequence
        is only a shallow rather than deep copy of the passed sequence.
    '''
    assert types.is_sequence_nonstr(omit_items), (
        types.assert_not_sequence_nonstr(omit_items))

    # Return a generator-based shallow copy of this sequence.
    return get_items_satisfying(
        sequence=sequence,
        item_satisfier=lambda item: item in omit_items
    )

# ....................{ REPLACERS                          }....................
def replace_items(
    sequence: "collections.Sequence",
    item_replacements: "collections.Mapping") -> (
    'collections.Sequence'):
    '''
    Get a new non-string sequence transformed from the passed non-string
    sequence by replacing all elements of the latter equalling any key of the
    passed dictionary by that key's value.

    This method effectively performs an equality-based global
    search-and-replacement for sequence elements.

    Parameters
    ----------
    sequence : collections.Sequence
        Original sequence to be returned transformed. For safety, this
        function does _not_ modify this sequence.
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
