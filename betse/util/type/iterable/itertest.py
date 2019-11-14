#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level **non-string iterable tester** (i.e., utility functions testing and
validating non-string objects implementing the abstract
:class:`collections.abc.Iterable` base class) facilities.
'''

# ....................{ IMPORTS                           }....................
from betse.exceptions import BetseIterableException
from betse.util.type.types import (
    type_check, CallableTypes, IterableTypes, SetType, TestableTypes)

# ....................{ EXCEPTIONS                        }....................
@type_check
def die_unless_items_instance_of(
    iterable: IterableTypes, cls: TestableTypes) -> None:
    '''
    Raise an exception unless *all* items of the passed iterable are instances
    of the passed class or tuple of classes.

    Parameters
    ----------
    iterable : IterableTypes
        Iterable to be validated.
    cls : TestableTypes
        Class or tuple of classes to validate that all items of this iterable
        be instances of.

    Raises
    ----------
    BetseIterableException
        If at least one item of this iterable is *not* an instance of this
        class or tuple of classes.

    See Also
    ----------
    :func:`is_items_instance_of`
        Further details.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.iterable import iterget

    # If one or more items of this iterable are *NOT* such instances...
    if not is_items_instance_of(iterable=iterable, cls=cls):
        # First such item.
        item_invalid = iterget.get_item_first_not_instance_of(
            iterable=iterable, cls=cls)

        # Raise an exception embedding this item.
        raise BetseIterableException(
            'Iterable item {} not instance of {}.'.format(item_invalid, cls))


@type_check
def die_unless_items_unique(iterable: IterableTypes) -> None:
    '''
    Raise an exception unless *all* items of the passed iterable are **unique**
    (i.e., no two distinct items are equal).

    Equivalently, this function raises an exception if one or more items of the
    passed iterable are duplicates and hence non-unique.

    Parameters
    ----------
    iterable : IterableTypes
        Iterable to be validated.

    Raises
    ----------
    BetseIterableException
        If one or more items of this iterable are duplicates.

    See Also
    ----------
    :func:`is_items_unique`
        Further details.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.iterable import generators, iterget
    from betse.util.type.text.string import strjoin

    # Iterable to be tested.
    #
    # If this iterable is a generator, coerce this generator into a tuple
    # *BEFORE* passing this generator to the is_items_unique() function, which
    # would internally consume and hence reduce this generator to the empty
    # generator, which would prevent the subsequent logic from detecting which
    # duplicate items of the original generator.
    iterable_testable = generators.to_tuple_if_generator(iterable)

    # If one or more items of this iterable are duplicates...
    if not is_items_unique(iterable_testable):
        # Set of all such duplicates.
        items_duplicate = iterget.get_items_duplicate(iterable_testable)

        # Human-readable conjunction of these duplicates, implicitly
        # converting each item of this set to a string as needed.
        items_duplicate_text = (
            strjoin.join_iterable_as_conjunction_double_quoted(
                items_duplicate))

        # Raise an exception embedding this set.
        raise BetseIterableException(
            'Iterable items non-unique: {}'.format(items_duplicate_text))

# ....................{ TESTERS                           }....................
@type_check
def is_reversible(iterable: IterableTypes) -> bool:
    '''
    ``True`` only if the passed iterable is **reversible** (i.e., successfully
    passable to the :func:`reversed` builtin).

    Specifically, this function returns ``True`` only if this iterable either:

    * Defines the ``__reversed__()`` special method.
    * Defines the ``__len__()`` *and* ``__getitem__()`` special methods, thus
      satisfying the :class:`collections.abc.Sequence` protocol.

    Parameters
    ----------
    iterable : IterableTypes
        Iterable to be inspected.

    Returns
    ----------
    bool
        ``True`` only if this iterable is reversible.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.obj import objtest

    # Return true only if this iterable defines either...
    return (
        # The __reversed__() special method *OR*...
        objtest.has_method(iterable, '__reversed__') or
        # The __len__() *AND* __getitem__() special methods.
        objtest.has_method(iterable, '__len__', '__getitem__'))

# ....................{ TESTERS ~ items                   }....................
@type_check
def is_items_instance_of(iterable: IterableTypes, cls: TestableTypes) -> bool:
    '''
    ``True`` only if *all* items of the passed iterable are instances of the
    passed class or tuple of classes.

    Parameters
    ----------
    iterable : IterableTypes
        Iterable to be searched.
    cls : TestableTypes
        Class or tuple of classes to test that all items of this iterable be
        instances of.

    Returns
    ----------
    bool
        ``True`` only if *all* items of this iterable are instances of this
        class or tuple of classes.
    '''

    return all(isinstance(item, cls) for item in iterable)


#FIXME: This function currently assumes all items of this iterable to be
#hashable (i.e., implement the "collections.abc.Hashable" interface). If this
#is *NOT* the case, exceptions and/or non-deterministic behaviour is likely.
#Ergo, this function must be generalized to:
#
#* Efficiently detect whether this iterable contains one or more non-hashable
#  items.
#* In that case, either:
#  * Raise an exception.
#  * Reduce to a less efficient approach supporting non-hashable items.
#    Notably, while "set" objects require hashability, "tuple" objects do not.
#    Ergo, leveraging a "tuple" rather than "set" should suffice here.
@type_check
def is_items_unique(iterable: IterableTypes) -> bool:
    '''
    ``True`` only if *all* items of the passed iterable are **unique** (i.e.,
    no two distinct items are equal).

    Specifically, this function:

    * If this iterable implements the :class:`collections.abc.Set` interface,
      returns ``True`` immediately.
    * Else if this iterable is sufficiently small (e.g., contains less than
      approximately 50 items), reduces this iterable to a :class:`set` and
      compares the length of these two containers for equality.
    * Else, efficiently iterates this iterable at most once and hence has
      worst-case time complexity ``O(len(iterable))``.

    Parameters
    ----------
    iterable : IterableTypes
        Iterable to be inspected.

    See Also
    ----------
    https://stackoverflow.com/a/5281641/2809027
        StackOverflow post strongly inspiring this implementation.

    Returns
    ----------
    bool
        ``True`` only if *all* items of this iterable are unique.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.obj import objtest

    # If this iterable is a set and hence guarantees uniqueness, efficiently
    # reduce to a noop.
    if isinstance(iterable, SetType):
        return True
    # Else, this iterable is *NOT* a set and hence does *NOT* guarantee
    # uniqueness.

    # If this iterable implements the "collections.abc.Sized" interface and
    # hence defines the __len__ dunder method implicitly called by the len()
    # builtin explicitly called here...
    if objtest.has_method(iterable, '__len__'):
        # Number of items in this iterable.
        iterable_len = len(iterable)

        # If this iterable is sufficiently small (where "sufficiently small"
        # has yet to be quantitatively measured and hence is qualitatively
        # arbitrary), return true only if this iterable contains exactly as
        # many items as the corresponding set of the same items, in which case
        # no such item is a duplicate of any other such item.
        #
        # Note that substituting "frozenset" for "set" here would yield no
        # tangible improvements to space or time complexity. The current
        # "frozenset" implementation is sadly naive and hence fails to exploit
        # obvious optimization opportunities (e.g., generation of a perfect
        # hash function specific to each "frozenset"). See also this pertinent
        # StackOverflow thread:
        #     https://stackoverflow.com/questions/36555214/set-vs-frozenset-performance
        if iterable_len < 64:
            return iterable_len == len(set(iterable))
        # Else, this iterable is sufficiently large to warrant an efficient
        # approach guaranteed to short circuit (i.e., "early exit") on the
        # first non-unique item of this iterable.

    # Set of all unique items of this iterable previously visited below.
    items_unique = set()

    # Return true only if no items in this iterable are duplicates (i.e., have
    # already been visited by a prior iteration of this test).
    return not any(
        # If this item is unique, add this item to this set as a side effect...
        item in items_unique or items_unique.add(item)
        # For each item of this iterable.
        for item in iterable)

# ....................{ TESTERS ~ item                    }....................
@type_check
def is_item_satisfying(
    iterable: IterableTypes, predicate: CallableTypes) -> bool:
    '''
    ``True`` only if some item of the passed iterable satisfies the passed
    **predicate** (i.e., callable accepting one parameter returning ``True``
    only if this parameter suffices).

    Parameters
    ----------
    iterable : IterableTypes
        Iterable to be searched.
    predicate : CallableTypes
        Callable accepting one parameter and returning ``True`` only if this
        parameter suffices.

    Returns
    ----------
    bool
        ``True`` only if some item of this iterable satisfies this predicate.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.iterable import iterget
    from betse.util.type.obj.sentinels import SENTINEL

    # First item satisfying this predicate in this iterable if any *OR* the
    # sentinel placeholder otherwise.
    first_item = iterget.get_item_first_satisfying_or_sentinel(
        iterable, predicate)

    # Return True only if this item is *NOT* the sentinel, in which case some
    # item satisfies this predicate.
    return first_item is not SENTINEL


@type_check
def is_item_instance_of(
    iterable: IterableTypes, cls: TestableTypes) -> bool:
    '''
    ``True`` only if some item of the passed iterable is an instance of the
    passed class or tuple of classes.

    Parameters
    ----------
    iterable : IterableTypes
        Iterable to be searched.
    cls : TestableTypes
        Class or tuple of classes of the item to search this iterable for.

    Returns
    ----------
    bool
        ``True`` only if some item of this iterable is an instance of this
        class or tuple of classes.
    '''

    return is_item_satisfying(
        iterable=iterable, predicate=lambda item: isinstance(item, cls))
