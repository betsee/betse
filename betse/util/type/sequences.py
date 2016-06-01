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
from collections.abc import Container, Mapping, Sequence

# ....................{ GETTERS                            }....................
def get_items_satisfying(
    sequence: Sequence, item_satisfier: callable) -> Sequence:
    '''
    New non-string sequence containing only the proper subset of elements from
    the passed non-string sequence satisfying the passed callable.

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
def get_items_prefixed_by(sequence: Sequence, item_prefix: str) -> Sequence:
    '''
    New non-string sequence containing only the proper subset of elements from
    the passed non-string sequence that are strings prefixed by the passed
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
def omit_item(sequence: Sequence, item: object) -> Sequence:
    '''
    New non-string sequence containing all elements of the first passed
    non-string sequence _not_ equalling the passed object.

    Parameters
    ----------
    sequence : collections.Sequence
        Original sequence to return a proper subset of. For safety, this
        function does _not_ modify this sequence.
    item : object
        Object to be omitted from the returned sequence.

    Returns
    ----------
    collections.Sequence
        Proper subset of the passed sequence. For efficiency, this new sequence
        is only a shallow rather than deep copy of the passed sequence.
    '''

    return omit_items(sequence=sequence, items=(item,))


def omit_items(sequence: Sequence, items: Container) -> Sequence:
    '''
    New non-string sequence containing all elements of the passed non-string
    sequence _not_ contained in the passed non-string container.

    Parameters
    ----------
    sequence : collections.Sequence
        Original sequence to return a proper subset of. For safety, this
        function does _not_ modify this sequence.
    items : collections.Container
        Container containing all elements to be omitted from the returned
        sequence.

    Returns
    ----------
    collections.Sequence
        Proper subset of the passed sequence. For efficiency, this new sequence
        is only a shallow rather than deep copy of the passed sequence.
    '''
    assert types.is_container_nonstr(items), (
        types.assert_not_container_nonstr(items))

    # Return a generator-based shallow copy of this sequence.
    return get_items_satisfying(
        sequence=sequence,
        item_satisfier=lambda item: item not in items,
    )

# ....................{ REMOVERS                           }....................
def remove_item(sequence: Sequence, item: object) -> None:
    '''
    Remove all elements equalling the passed object from the passed non-string
    sequence.

    Parameters
    ----------
    sequence : collections.Sequence
        Sequence to remove elements from.
    item : object
        Object to be removed.
    '''

    remove_items(sequence=sequence, items=(item,))


def remove_items(sequence: Sequence, items: Container) -> None:
    '''
    Remove all elements contained in the passed non-string container from
    the passed non-string sequence.

    Parameters
    ----------
    sequence : collections.Sequence
        Sequence to remove elements from.
    items : collections.Container
        Container containing all elements to be remove from this sequence.
    '''
    assert types.is_sequence_nonstr(sequence), (
        types.assert_not_sequence_nonstr(sequence))

    # Slice assignment implicitly modifies the original sequence, permitting
    # efficient reuse of existing assignment-based functionality.
    sequence[:] = omit_items(sequence, items)

# ....................{ REPLACERS                          }....................
#FIXME: Rename the "item_replacements" parameter to merely "replacements".
def replace_items(sequence: Sequence, item_replacements: Mapping) -> Sequence:
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
def sort_lexicographic_ascending(sequence: Sequence) -> Sequence:
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
