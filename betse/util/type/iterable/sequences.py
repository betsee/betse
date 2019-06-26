#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level **non-string sequence** (i.e., non-string object implementing the
abstract base class `collections.abc.Sequence`) facilities.

See Also
----------
:func:`betse.util.type.types.is_sequence`
    Further details on what constitutes sequences and non-string sequences.
'''

#FIXME: Most if not all of the existing functionality provided by this
#submodule should be generalized to support higher-level iterable rather than
#lower-level sequence types and then shifted into the existing "iterables"
#submodule.

# ....................{ IMPORTS                           }....................
from betse.exceptions import BetseSequenceException
from betse.util.type import types
from betse.util.type.types import (
    type_check,
    CallableTypes,
    ContainerTypes,
    IterableTypes,
    MappingType,
    SequenceTypes,
    StrOrNoneTypes,
)

# ....................{ EXCEPTIONS                        }....................
@type_check
def die_unless_index(sequence: SequenceTypes, index: int) -> None:
    '''
    Raise an exception unless the passed index **indexes** (i.e., is a
    non-negative integer strictly greater than -1 and less than the length of)
    the passed sequence.

    Equivalently, this function raises an exception this index is either:

    * Less than 0.
    * Greater than or equal to the length of this sequence.

    Parameters
    ----------
    sequence : SequenceTypes
        Sequence to be validated.
    index : int
        0-based index to validate against this sequence.

    Raises
    ----------
    BetseSequenceException
        Unless this index indexes this sequence.
    '''

    # If this index does *NOT* index this sequence, raise an exception.
    if not is_index(sequence=sequence, index=index):
        raise BetseSequenceException(
            'Sequence index {} invalid (i.e., not in range [0, {}]).'.format(
                index, len(sequence)))

# ....................{ EXCEPTIONS ~ length               }....................
@type_check
def die_if_empty(
    sequence: SequenceTypes, exception_message: StrOrNoneTypes = None) -> None:
    '''
    Raise an exception with the passed message (defaulting to a human-readable
    message) if the passed sequence is **empty** (i.e., contains *no* items).

    Parameters
    ----------
    sequence: SequenceTypes
        Sequence to be validated.
    exception_message : StrOrNoneTypes
        Exception message to be raised. Defaults to ``None``, in which case an
        exception message synthesized from the passed arguments is raised.

    Raises
    ----------
    BetseSequenceException
        If the sequence is empty.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.obj import objects

    # If this sequence is empty, raise an exception.
    if is_empty(sequence):
        # If no exception message was passed, synthesize one.
        if not exception_message:
            exception_message = 'Sequence "{}" empty.'.format(
                objects.get_class_name_unqualified(sequence))

        # Raise this exception.
        raise BetseSequenceException(exception_message)



@type_check
def die_if_length_less_than(
    # Mandatory parameters.
    sequence: SequenceTypes,
    length: int,

    # Optional parameters.
    exception_message: StrOrNoneTypes = None,
) -> None:
    '''
    Raise an exception with the passed message (defaulting to a human-readable
    message) unless the length of the passed sequence is greater than or equal
    to the passed length (i.e., unless this sequence contains at least as many
    items as this length).

    Equivalently, this function raises an exception if this sequence contains
    strictly less than this many items.

    Parameters
    ----------
    sequence: SequenceType
        Sequence to be validated.
    length: int
        Minimum sequence length to test for.
    exception_message : StrOrNoneTypes
        Exception message to be raised. Defaults to ``None``, in which case an
        exception message synthesized from the passed arguments is raised.

    Raises
    ----------
    BetseSequenceException
        If this sequence contains strictly less than this many items.
    '''

    # If this sequence contains less than this many items, raise an exception.
    if len(sequence) < length:
        # If no exception message was passed, synthesize one.
        if not exception_message:
            exception_message = 'Sequence length {} < {}.'.format(
                len(sequence), length)

        # Raise this exception.
        raise BetseSequenceException(exception_message)


@type_check
def die_unless_length(
    # Mandatory parameters.
    sequence: SequenceTypes,
    length: int,

    # Optional parameters.
    exception_message: StrOrNoneTypes = None,
) -> None:
    '''
    Raise an exception with the passed message (defaulting to a human-readable
    message) unless the passed sequence is of the passed length.

    Equivalently, this function raises an exception if this sequence is *not*
    of this length.

    Parameters
    ----------
    sequence: SequenceType
        Sequence to be validated.
    length: int
        Exact sequence length to test for.
    exception_message : StrOrNoneTypes
        Exception message to be raised. Defaults to ``None``, in which case an
        exception message synthesized from the passed arguments is raised.

    Raises
    ----------
    BetseSequenceException
        If this sequence is *not* of this length.
    '''

    # If this sequence is *NOT* of this length, raise an exception.
    if len(sequence) != length:
        # If no exception message was passed, synthesize one.
        if not exception_message:
            exception_message = 'Sequence length {} not {}.'.format(
                len(sequence), length)

        # Raise this exception.
        raise BetseSequenceException(exception_message)

# ....................{ TESTERS                           }....................
@type_check
def is_empty(*sequences: SequenceTypes) -> bool:
    '''
    ``True`` only if all passed sequences are **empty** (i.e., contain no
    elements).

    Parameters
    ----------
    sequences : tuple[SequenceTypes]
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


@type_check
def is_index(sequence: SequenceTypes, index: int) -> bool:
    '''
    ``True`` only if the passed index **indexes** (i.e., is a non-negative
    integer strictly greater than -1 and less than the length of) the passed
    sequence.

    Parameters
    ----------
    sequence : SequenceTypes
        Sequence to be tested.
    index : int
        0-based index to test against this sequence.

    Returns
    ----------
    bool
        ``True`` only if this index indexes this sequence.
    '''

    return 0 <= index < len(sequence)

# ....................{ TESTERS ~ type                    }....................
def is_sequence(*objs: object) -> bool:
    '''
    ``True`` only if all passed objects are **sequences** (i.e., of types
    conforming to but *not* necessarily subclassing the canonical
    :class:`collections.abc.Sequence` API).

    Parameters
    ----------
    objs : tuple[object]
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


#FIXME: Excise everywhere in favour of the more appropriate
#betse.lib.numpy.nparray.is_array() function -- which, admittedly, may need to
#be generalized to accept variadic parameters.
def is_numpy_array(*objs: object) -> bool:
    '''
    ``True`` only if all passed objects are **Numpy arrays** (i.e., instances
    of the :class:`numpy.ndarray` superclass).

    Parameters
    ----------
    objs : tuple[object]
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

# ....................{ GETTERS                           }....................
@type_check
def get_index(sequence: SequenceTypes, index: int) -> object:
    '''
    Item with the passed index in this sequence if this index **indexes**
    (i.e., is a non-negative integer strictly greater than -1 and less than the
    length of) this sequence *or*  raise an exception otherwise (i.e., if this
    index does *not* index this sequence).

    Parameters
    ----------
    sequence : SequenceTypes
        Sequence to be inspected.
    index : int
        0-based index to return the item in this sequence of.

    Returns
    ----------
    object
        Item with this index in this sequence.

    Raises
    ----------
    BetseSequenceException
        If this index does *not* index this sequence.
    '''

    # If this index does *NOT* index this sequence, raise an exception.
    die_unless_index(sequence=sequence, index=index)
    # Else, this index indexes this sequence.

    # Return the item with this index in this sequence.
    return sequence[index]

# ....................{ GETTERS ~ items                   }....................
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

        * ``True`` if this element satisfies the desired requirements.
        * ``False`` otherwise.

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
        function does *not* modify this sequence.
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
            types.is_str(item) and item.startswith(item_prefix)))

# ....................{ OMITTERS                          }....................
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
def omit_items(
    sequence: SequenceTypes, items: ContainerTypes) -> SequenceTypes:
    '''
    New non-string sequence containing all items of the passed non-string
    sequence *not* contained in the passed non-string container.

    Parameters
    ----------
    sequence : SequenceTypes
        Original sequence to return a proper subset of. For safety, this
        function does *not* modify this sequence.
    items : ContainerTypes
        Container of all items to be omitted from the returned sequence.

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

# ....................{ REMOVERS                          }....................
def remove_item(sequence: SequenceTypes, item: object) -> None:
    '''
    Remove all items equal to the passed object from the passed mutable
    non-string sequence in-place.

    Parameters
    ----------
    sequence : SequenceTypes
        Sequence to remove items from in-place.
    item : object
        Object to be removed.
    '''

    remove_items(sequence=sequence, items=(item,))


@type_check
def remove_items(sequence: SequenceTypes, items: ContainerTypes) -> None:
    '''
    Remove all items contained in the passed non-string container from the
    passed mutable non-string sequence in-place.

    Parameters
    ----------
    sequence : SequenceTypes
        Sequence to remove items from in-place.
    items : ContainerTypes
        Container of all items to be removed from this sequence.
    '''

    # Slice assignment implicitly modifies the original sequence, permitting
    # efficient reuse of existing assignment-based functionality.
    sequence[:] = omit_items(sequence, items)

# ....................{ REPLACERS                         }....................
@type_check
def replace_items(
    sequence: SequenceTypes, replacements: MappingType) -> SequenceTypes:
    '''
    New non-string sequence containing all items of the passed non-string
    sequence (_in the same order_) such that each item equal to a key of the
    passed dictionary is replaced by the value of that key.

    This method effectively performs an equality-based global
    search-and-replacement for sequence items.

    Parameters
    ----------
    sequence : SequenceTypes
        Original sequence to be returned transformed. For safety, this
        function does *not* modify this sequence.
    replacements : MappingType
        Mapping whose:

        * Keys are input values to find in the passed sequence. For simplicity,
          keys are matched via object equality rather than more complex object
          matching alternatives (e.g., glob, regex, prefix, suffix, substring).
        * Values are output values to replace these input values with in the
          returned sequence.

    Returns
    ----------
    SequenceTypes
        Sequence transformed from the passed sequence. For efficiency, this
        sequence is only a shallow rather than deep copy of the passed
        sequence.
    '''

    # Type of both the passed sequence and the sequence to be returned.
    sequence_type = type(sequence)

    # Return a shallow copy of this sequence, replacing each item that is a key
    # of the passed mapping by that key's value.
    return sequence_type(
        replacements[item] if item in replacements else item
        for item in sequence
    )

# ....................{ CONVERTERS                        }....................
@type_check
def to_sequence(iterable: IterableTypes) -> SequenceTypes:
    '''
    Sequence converted from the passed iterable if this iterable is not a
    sequence *or* this iterable as is otherwise (i.e., if this iterable is a
    sequence).

    For efficiency, this function only performs a shallow copy of the items in
    this iterable on creating this sequence. Ergo, this sequence is guaranteed
    to contain the same items in the same order as this iterable.

    Parameters
    ----------
    iterable : IterableTypes
        Iterable to be converted into a sequence.

    Returns
    ----------
    SequenceTypes
        Either:

        * If this iterable is already a sequence, this sequence as is.
        * Else, a new sequence converted from this iterable (and hence
          containing the same items in the same order as this iterable).
    '''

    # One-liners to rule them all.
    return iterable if is_sequence(iterable) else tuple(iterable)
