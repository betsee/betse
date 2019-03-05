#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level **non-string iterable testers** (i.e., utility functions testing and
validating non-string objects implementing the abstract
:class:`collections.abc.Iterable` base class) facilities.
'''

# ....................{ IMPORTS                           }....................
from betse.exceptions import BetseIterableException
from betse.util.type.types import (
    type_check, CallableTypes, IterableTypes, TestableTypes)

# ....................{ EXCEPTIONS                        }....................
@type_check
def die_unless_items_instance_of(
    iterable: IterableTypes, cls: TestableTypes) -> None:
    '''
    Raise an exception unless *all* items of the passed iterable are instances
    of the passed class or tuple of classes.

    Parameters
    ----------
    iterable: IterableTypes
        Iterable to be validated.
    cls : TestableTypes
        Class or tuple of classes to validate that all items of this iterable
        be instances of.

    Raises
    ----------
    BetseIterableException
        If at least one item of this iterable is *not* an instance of this
        class or tuple of classes.
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

    Parameters
    ----------
    iterable: IterableTypes
        Iterable to be validated.

    Raises
    ----------
    BetseIterableException
        If at least one item of this iterable is a duplicate.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.iterable import iterget

    # If one or more items of this iterable are duplicates...
    if not is_items_unique(iterable):
        # Set of all such duplicates.
        items_duplicate = iterget.get_items_duplicate(iterable)

        # Raise an exception embedding this set.
        raise BetseIterableException(
            'Iterable items {} duplicate.'.format(items_duplicate))

# ....................{ TESTERS                           }....................
@type_check
def is_reversible(iterable: IterableTypes) -> bool:
    '''
    ``True`` only if the passed iterable is **reversible** (i.e., successfully
    passable to the :func:`reversed` builtin).

    Specifically, this function returns ``True`` only if this iterable either:

    * Defines the ``__reversed__()`` special method.
    * Defines the ``__len__()`` *and* ``__getitem__()`` special methods,
      satisfying the sequence protocol.

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
    from betse.util.type.obj import objects

    # Return True only if this iterable defines either...
    return (
        # The __reversed__() special method *OR*...
        objects.has_method(iterable, '__reversed__') or
        # The __len__() *AND* __getitem__() special methods.
        objects.has_method(iterable, '__len__', '__getitem__'))

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


@type_check
def is_items_unique(iterable: IterableTypes) -> bool:
    '''
    ``True`` only if *all* items of the passed iterable are **unique** (i.e.,
    no two distinct items are equal).

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

    # Set of all unique items of this iterable previously visited below.
    items_unique = set()

    # Return True only if no items in this iterable are duplicates (i.e., have
    # already been visited by a prior iteration of this test).
    return not any(
        # If this item is unique, add this item to this set as a side effect.
        item in items_unique or items_unique.add(item)
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
