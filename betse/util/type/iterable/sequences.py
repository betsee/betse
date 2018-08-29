#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level **non-string sequence** (i.e., non-string object implementing the
abstract base class `collections.abc.Sequence`) facilities.

See Also
----------
betse.util.type.types.is_sequence
    Further details on what constitutes sequences and non-string sequences.
'''

#FIXME: Most if not all of the existing functionality provided by this submodule
#should be generalized to support higher-level iterable rather than lower-level
#sequence types and then shifted into the existing "iterables" submodule.
#Ideally, this submodule could then be excised entirely.

# ....................{ IMPORTS                            }....................
from betse.exceptions import BetseSequenceException
from betse.util.type import types
from betse.util.type.types import (
    type_check, CallableTypes, SequenceTypes, StrOrNoneTypes)
from collections.abc import Container, Mapping

# ....................{ EXCEPTIONS                         }....................
@type_check
def die_if_empty(
    *sequences: SequenceTypes, label: str = 'Sequence') -> None:
    '''
    Raise an exception prefixed by the passed label unless all passed sequences
    are **non-empty** (i.e., contain at least one element).

    Parameters
    ----------
    sequences: tuple
        Tuple of all sequences to be validated.
    label : optional[str]
        Human-readable label prefixing exception messages raised by this method.
        Defaults to a general-purpose string.

    Raises
    ----------
    BetseSequenceException
        If any passed sequence is empty.
    '''

    # If only one sequence is passed...
    if len(sequences) == 1:
        # If this sequence is non-empty, raise a simplistic exception.
        if is_empty(sequences[0]):
            raise BetseSequenceException('{} empty.'.format(label.capitalize()))
    # Else, multiple sequences are passed.
    else:
        # For each such sequence...
        for sequence_index, sequence in enumerate(sequences):
            # If this sequence is non-empty, raise an exception identifying the
            # index of this sequence in this parameter list.
            if is_empty(sequence):
                raise BetseSequenceException(
                    '{} {} empty.'.format(label.capitalize(), sequence_index))


@type_check
def die_unless_len(
    sequence: SequenceTypes,
    sequence_len: int,
    exception_message: StrOrNoneTypes = None,
) -> None:
    '''
    Raise an exception with the passed message (defaulting to a human-readable
    message) if the passed sequence is *not* of the passed length.

    Parameters
    ----------
    sequence: SequenceType
        Sequence to be validated.
    len: int
        Sequence length to test for.
    exception_message : optional[str]
        Exception message to be raised. Defaults to ``None``, in which case an
        exception message synthesized from the passed arguments is raised.

    Raises
    ----------
    BetseSequenceException
        If this sequence is *not* of this length.
    '''

    # If this sequence is *NOT* of this length, raise an exception.
    if len(sequence) != sequence_len:
        # If no exception message was passed, synthesize one.
        if not exception_message:
            exception_message = 'Sequence length {} not {}.'.format(
                len(sequence), sequence_len)

        # Raise this exception.
        raise BetseSequenceException(exception_message)

# ....................{ TESTERS ~ type                     }....................
def is_sequence(*objs: object) -> bool:
    '''
    ``True`` only if all passed objects are **sequences** (i.e., of types
    conforming to but *not* necessarily subclassing the canonical
    :class:`collections.abc.Sequence` API).

    Parameters
    ----------
    objs: tuple[object]
        Tuple of all objects to be tested.

    Returns
    ----------
    bool
        ``True`` only if these objects are all sequences.

    See Also
    ----------
    :class:`collections.abc.Sequence`
        Further details.
    '''

    # all(). It is awesome.
    return all(isinstance(obj, SequenceTypes) for obj in objs)


def is_numpy_array(*objs: object) -> bool:
    '''
    ``True`` only if all passed objects are **Numpy arrays** (i.e., instances of
    the :class:`numpy.ndarray` superclass).

    Parameters
    ----------
    objs: tuple[object]
        Tuple of all objects to be tested.

    Returns
    ----------
    bool
        ``True`` only if these objects are all Numpy arrays.
    '''

    # Avoid importing third-party packages at the top level, for safety.
    from numpy import ndarray

    # all(). It exceeds at all the tests.
    return all(isinstance(obj, ndarray) for obj in objs)

# ....................{ TESTERS ~ len                      }....................
@type_check
def is_empty(*sequences: SequenceTypes) -> bool:
    '''
    ``True`` only if all passed sequences are **empty** (i.e., contain no
    elements).

    Parameters
    ----------
    sequences: tuple
        Tuple of all sequences to be tested.

    Returns
    ----------
    bool
        ``True`` only if these sequences are all empty.
    '''

    # all(). It is all awesome.
    #
    # To transparently support Numpy arrays, the length of this sequence *MUST*
    # be explicitly tested via the len() builtin rather than implicitly tested
    # as a boolean. On the latter, Numpy raises the following exception:
    #
    #     ValueError: The truth value of an array with more than one element is
    #     ambiguous. Use a.any() or a.all()
    return all(len(sequence) == 0 for sequence in sequences)

# ....................{ GETTERS                            }....................
#FIXME: Generalize into a function of the same name accepting an iterable rather
#than sequence and shift into the "iterables" submodule. When doing so, rename
#the "item_satisfier" parameter to "predicate" for orthogonality and clarity.
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
        CallableTypes (e.g., function, lambda) accepting a single element of
        this sequence and returning only:
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
#FIXME: Generalize into a function of the same name accepting an iterable rather
#than sequence and shift into the "iterables" submodule. When doing so, rename
#the "item_prefix" parameter to simply "prefix".
@type_check
def get_items_prefixed_by(
    sequence: SequenceTypes, item_prefix: str) -> SequenceTypes:
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
    non-string sequence *not* equalling the passed object.

    Parameters
    ----------
    sequence : SequenceTypes
        Original sequence to return a proper subset of. For safety, this
        function does *not* modify this sequence.
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
    sequence *not* contained in the passed non-string container.

    Parameters
    ----------
    sequence : SequenceTypes
        Original sequence to return a proper subset of. For safety, this
        function does *not* modify this sequence.
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
    Remove all elements equal to the passed object from the passed mutable
    non-string sequence in-place.

    Parameters
    ----------
    sequence : collections.Sequence
        Sequence to remove elements from in-place.
    item : object
        Object to be removed.
    '''

    remove_items(sequence=sequence, items=(item,))


@type_check
def remove_items(sequence: SequenceTypes, items: Container) -> None:
    '''
    Remove all elements contained in the passed non-string container from
    the passed mutable non-string sequence in-place.

    Parameters
    ----------
    sequence : SequenceTypes
        Sequence to remove elements from in-place.
    items : Container
        Container containing all elements to be remove from this sequence.
    '''

    # Slice assignment implicitly modifies the original sequence, permitting
    # efficient reuse of existing assignment-based functionality.
    sequence[:] = omit_items(sequence, items)

# ....................{ REPLACERS                          }....................
@type_check
def replace_items(
    sequence: SequenceTypes, replacements: Mapping) -> SequenceTypes:
    '''
    New non-string sequence containing all elements of the passed non-string
    sequence (_in the same order_) such that each element equal to a key of the
    passed dictionary is replaced by the value of that key.

    This method effectively performs an equality-based global
    search-and-replacement for sequence elements.

    Parameters
    ----------
    sequence : collections.Sequence
        Original sequence to be returned transformed. For safety, this
        function does _not_ modify this sequence.
    replacements : collections.Mapping
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

    # Type of both the passed sequence and the sequence to be returned.
    sequence_type = type(sequence)

    # Return a shallow copy of this sequence, replacing each element that is a
    # key of the passed mapping by that key's value.
    return sequence_type(
        replacements[item] if item in replacements else item
        for item in sequence
    )
