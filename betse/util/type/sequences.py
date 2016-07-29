#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level **non-string sequence** (i.e., non-string object implementing the
abstract base class `collections.abc.Sequence`) facilities.

See Also
----------
betse.util.type.types.is_sequence
    Further details on what constitutes sequences and non-string sequences.
'''

# ....................{ IMPORTS                            }....................
from betse.util.type import types
from betse.util.type.types import type_check, CallableTypes, SequenceTypes
from collections.abc import Container, Mapping

# ....................{ GETTERS                            }....................
@type_check
def get_items_satisfying(
    sequence: SequenceTypes, item_satisfier: CallableTypes) -> SequenceTypes:
    '''
    New non-string sequence containing only the proper subset of elements from
    the passed non-string sequence satisfying the passed callable.

    This method effectively performs a general-purpose global search for
    sequence elements.

    Parameters
    ----------
    sequence : SequenceTypes
        Original sequence to return a proper subset of. For safety, this
        function does _not_ modify this sequence.
    item_satisfier : CallableTypes
        CallableTypes (e.g., function, lambda) accepting a single element of this
        sequence and returning only:
        * `True` if this element satisfies the desired requirements.
        * `False` otherwise.

    Returns
    ----------
    collections.Sequence
        Proper subset of the passed sequence. For efficiency, this new sequence
        is only a shallow rather than deep copy of the passed sequence.
    '''

    # Type of both the passed sequence and the sequence to be returned.
    sequence_type = type(sequence)

    # Return a generator-based shallow copy of this sequence.
    return sequence_type(item for item in sequence if item_satisfier(item))

# ....................{ GETTERS ~ str                      }....................
@type_check
def get_items_prefixed_by(sequence: SequenceTypes, item_prefix: str) -> SequenceTypes:
    '''
    New non-string sequence containing only the proper subset of elements from
    the passed non-string sequence that are strings prefixed by the passed
    string prefix.

    This method effectively performs a prefix-based global search for sequence
    elements.

    Parameters
    ----------
    sequence : SequenceTypes
        Original sequence to return a proper subset of. For safety, this
        function does _not_ modify this sequence.
    item_prefix : str
        String prefixing all elements of the returned sequence.

    Returns
    ----------
    SequenceTypes
        Proper subset of the passed sequence. For efficiency, this new sequence
        is only a shallow rather than deep copy of the passed sequence.
    '''

    # Return a generator-based shallow copy of this sequence.
    return get_items_satisfying(
        sequence=sequence,
        item_satisfier=lambda item: (
            types.is_str(item) and item.startswith(item_prefix)
        ),
    )

# ....................{ OMITTERS                           }....................
def omit_item(sequence: SequenceTypes, item: object) -> SequenceTypes:
    '''
    New non-string sequence containing all elements of the first passed
    non-string sequence _not_ equalling the passed object.

    Parameters
    ----------
    sequence : SequenceTypes
        Original sequence to return a proper subset of. For safety, this
        function does _not_ modify this sequence.
    item : object
        Object to be omitted from the returned sequence.

    Returns
    ----------
    SequenceTypes
        Proper subset of the passed sequence. For efficiency, this new sequence
        is only a shallow rather than deep copy of the passed sequence.
    '''

    return omit_items(sequence=sequence, items=(item,))


@type_check
def omit_items(sequence: SequenceTypes, items: Container) -> SequenceTypes:
    '''
    New non-string sequence containing all elements of the passed non-string
    sequence _not_ contained in the passed non-string container.

    Parameters
    ----------
    sequence : SequenceTypes
        Original sequence to return a proper subset of. For safety, this
        function does _not_ modify this sequence.
    items : Container
        Container containing all elements to be omitted from the returned
        sequence.

    Returns
    ----------
    SequenceTypes
        Proper subset of the passed sequence. For efficiency, this new sequence
        is only a shallow rather than deep copy of the passed sequence.
    '''

    # Return a generator-based shallow copy of this sequence.
    return get_items_satisfying(
        sequence=sequence,
        item_satisfier=lambda item: item not in items,
    )

# ....................{ REMOVERS                           }....................
def remove_item(sequence: SequenceTypes, item: object) -> None:
    '''
    Remove all elements equalling the passed object from the passed non-string
    sequence.

    Parameters
    ----------
    sequence : collections.Sequence
        SequenceTypes to remove elements from.
    item : object
        Object to be removed.
    '''

    remove_items(sequence=sequence, items=(item,))


@type_check
def remove_items(sequence: SequenceTypes, items: Container) -> None:
    '''
    Remove all elements contained in the passed non-string container from
    the passed non-string sequence.

    Parameters
    ----------
    sequence : SequenceTypes
        SequenceTypes to remove elements from.
    items : Container
        Container containing all elements to be remove from this sequence.
    '''

    # Slice assignment implicitly modifies the original sequence, permitting
    # efficient reuse of existing assignment-based functionality.
    sequence[:] = omit_items(sequence, items)

# ....................{ REPLACERS                          }....................
#FIXME: Rename the "item_replacements" parameter to merely "replacements".
@type_check
def replace_items(sequence: SequenceTypes, item_replacements: Mapping) -> SequenceTypes:
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
        SequenceTypes transformed from the passed sequence. For efficiency, this
        sequence is only a shallow rather than deep copy of the passed sequence.
    '''

    # Type of both the passed sequence and the sequence to be returned.
    sequence_type = type(sequence)

    # Return a shallow copy of this sequence, replacing each element that is a
    # key of the passed mapping by that key's value.
    return sequence_type(
        item_replacements[item] if item in item_replacements else item
        for item in sequence
    )

# ....................{ SORTERS                            }....................
@type_check
def sort_lexicographic_ascending(sequence: SequenceTypes) -> SequenceTypes:
    '''
    Get a new non-string sequence sorted from the passed non-string sequence in
    **ascending lexicographic order** (i.e., traditional order of dead-tree
    dictionaries and encyclopedias).

    Parameters
    ----------
    sequence : SequenceTypes
        Unsorted sequence to be returned sorted. For generality, this sequence
        is _not_ modified by this function.

    Returns
    ----------
    SequenceTypes
        SequenceTypes sorted from the passed sequence. For efficiency, this sequence
        is only a shallow rather than deep copy of the passed sequence.
    '''

    return sorted(sequence)   # Well, that was easy.
