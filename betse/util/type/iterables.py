#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level **non-string iterable** (i.e., non-string object implementing the
abstract base class `collections.abc.Iterable`) facilities.

See Also
----------
:func:`betse.util.type.types.is_iterable`
    Further details on what constitutes iterables and non-string iterables.
'''

# ....................{ IMPORTS                            }....................
import itertools
from betse.exceptions import BetseIterableException
from betse.util.type import types
from betse.util.type.types import (
    type_check,
    CallableTypes,
    ClassType,
    GeneratorType,
    IterableTypes,
    MappingType,
    SizedType,
    TestableTypes,
)
from collections import deque
from operator import itemgetter

# ....................{ EXCEPTIONS                         }....................
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
        Class or tuple of classes to validate that all items of this iterable be
        instances of.

    Raises
    ----------
    BetseIterableException
        If at least one item of this iterable is *not* an instance of this class
        or tuple of classes.
    '''

    # If one or more items of this iterable are *NOT* such instances...
    if not is_items_instance_of(iterable=iterable, cls=cls):
        # First such item.
        item_invalid = get_item_first_not_instance_of(
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

    # If one or more items of this iterable are duplicates...
    if not is_items_unique(iterable):
        # Set of all such duplicates.
        items_duplicate = get_items_duplicate(iterable)

        # Raise an exception embedding this set.
        raise BetseIterableException(
            'Iterable items {} duplicate.'.format(items_duplicate))

# ....................{ TESTERS                            }....................
@type_check
def is_reversible(iterable: IterableTypes) -> bool:
    '''
    `True` only if the passed iterable is **reversible** (i.e., successfully
    passable to the :func:`reversed` builtin).

    Specifically, this function returns `True` if this iterable either:

    * Defines the `__reversed__()` special method.
    * Defines the `__len__()` and `__getitem__()` special methods, satisfying
      the sequence protocol.

    Parameters
    ----------
    iterable : IterableTypes
        Iterable to be inspected.

    Returns
    ----------
    bool
        `True` only if this iterable is reversible.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.obj import objects

    # Return True only if this iterable either...
    return (
        # Defines the __reversed__() special method.
        objects.is_method(iterable, '__reversed__') or (
            # Defines the __len__() and __getitem__() special methods.
            objects.is_method(iterable, '__len__') and
            objects.is_method(iterable, '__getitem__')
        )
    )

# ....................{ TESTERS ~ items                    }....................
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

# ....................{ TESTERS ~ item                     }....................
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
    from betse.util.type.obj.sentinels import SENTINEL

    # First item satisfying this predicate in this iterable if any *OR* the
    # sentinel placeholder otherwise.
    first_item = get_item_first_satisfying_or_sentinel(iterable, predicate)

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

# ....................{ GETTERS                            }....................
@type_check
def get_items_duplicate(iterable: IterableTypes) -> set:
    '''
    Set of all duplicate items in the passed iterable.

    Parameters
    ----------
    iterable : IterableTypes
        Iterable to be inspected.

    Returns
    ----------
    set
        Set of all duplicate items in this iterable.

    See Also
    ----------
    https://stackoverflow.com/a/9835819/2809027
        Stackoverflow answer strongly inspiring this implementation.
    '''

    # Set of all unique items of this iterable previously visited below.
    items_unique = set()

    # Set of all duplicate items of this iterable to be returned.
    items_duplicate = set()

    # For each item of this iterable...
    for item in iterable:
        # If this item has been previously visited, this item is a duplicate.
        if item in items_unique:
            items_duplicate.add(item)
        # Else, this item is unique.
        else:
            items_unique.add(item)

    # Return this set of all duplicate items.
    return items_duplicate

# ....................{ GETTERS ~ first                    }....................
@type_check
def get_item_first(iterable: IterableTypes) -> object:
    '''
    First item non-destructively retrieved from the passed iterable if this
    iterable is non-empty *or* raise an exception otherwise (i.e., if this
    iterable is empty).

    If the passed iterable is a:

    * Sequence, this is guaranteed to be the first item of this sequence.
    * Non-sequence (e.g., :class:`set`, :class:`dict`), this should be assumed
      to be a random item. While most non-sequences guarantee predictable order
      of retrieval assuming no intervening changes, this is a fairly unreliable
      assumption.

    Parameters
    ----------
    iterable : IterableTypes
        Iterable to be inspected.

    Returns
    ----------
    object
        First element non-destructively retrieved from this iterable.

    Raises
    ----------
    BetseIterableException
        If this iterable is empty.

    See Also
    ----------
    https://stackoverflow.com/a/40054478/2809027
        Cecil's Stackoverflow answer strongly inspiring this implementation,
        complete with detailed timings of all alternative solutions.
    '''

    # If this iterable is empty, raise an exception. Since iteration of empty
    # iterables always succeeds, this condition must be manually tested
    # beforehand. Failing to do so would result Python raising the following
    # non-human-readable exception below:
    #
    #     NameError: name 'first_item' is not defined
    if not iterable:
        raise BetseIterableException('Iterable "{}" empty.'.format(iterable))

    # Yup! Shockingly, the most verbose and unwieldy solution is the fastest.
    # Break immediately after the first iteration of this iterable.
    first_item = None
    for first_item in iterable:
        break

    # Return the first element iterated above.
    return first_item

# ....................{ GETTERS ~ first : instance         }....................
@type_check
def get_item_first_instance_of(
    iterable: IterableTypes, cls: TestableTypes, **kwargs) -> object:
    '''
    First instance of the passed class or tuple of classes retrieved from the
    passed iterable if this iterable contains such an item *or* raise an exception
    otherwise (i.e., if this iterable contains no such item).

    Parameters
    ----------
    iterable : IterableTypes
        Iterable to be searched.
    cls : TestableTypes
        Class or tuple of classes of the item to search this iterable for.
    kwargs : dict
        Dictionary of all remaining keyword arguments to be passed as is to the
        :func:`get_item_first_satisfying` function.

    Returns
    ----------
    object
        First instance of this class in this iterable.

    Raises
    ----------
    BetseIterableException
        If this iterable contains no such item.

    See Also
    ----------
    :func:`get_item_first_satisfying_or_sentinel`
        Further details on ordering guarantees.
    '''

    return get_item_first_satisfying(
        iterable=iterable,
        predicate=lambda item: isinstance(item, cls),
        **kwargs
    )


@type_check
def get_item_first_not_instance_of(
    iterable: IterableTypes, cls: TestableTypes, **kwargs) -> object:
    '''
    First item of this iterable that is *not* an instance of the passed class or
    tuple of classes if this iterable contains such an item *or* raise an
    exception otherwise (i.e., if all items of this iterable are such
    instances).

    Parameters
    ----------
    iterable : IterableTypes
        Iterable to be searched.
    cls : TestableTypes
        Class or tuple of classes of the item to search this iterable for.
    kwargs : dict
        Dictionary of all remaining keyword arguments to be passed as is to the
        :func:`get_item_first_satisfying` function.

    Returns
    ----------
    object
        First non-instance of this class in this iterable.

    Raises
    ----------
    BetseIterableException
        If this iterable contains only such instances.

    See Also
    ----------
    :func:`get_item_first_satisfying_or_sentinel`
        Further details on ordering guarantees.
    '''

    return get_item_first_satisfying(
        iterable=iterable,
        predicate=lambda item: not isinstance(item, cls),
        **kwargs
    )

# ....................{ GETTERS ~ first : satisfying       }....................
@type_check
def get_item_first_satisfying_or_sentinel(
    iterable: IterableTypes, predicate: CallableTypes) -> object:
    '''
    First item of the passed iterable satisfying the passed **predicate** (i.e.,
    callable accepting one parameter, returning `True` only if this parameter
    suffices) if this iterable contains such an item *or* the
    :attr:`betse.util.type.obj.sentinels.SENTINEL` placeholder constant
    otherwise.

    If the passed iterable is a:

    * Sequence, this is guaranteed to be the first such element.
    * Non-sequence (e.g., :class:`set`, :class:`dict`), this should be assumed
      to be a random such element. While most non-sequences guarantee
      predictable order of retrieval assuming no intervening changes, this is a
      fairly unreliable assumption.

    Parameters
    ----------
    iterable : IterableTypes
        Iterable to be inspected.
    predicate : CallableTypes
        Callable accepting one parameter and returning ``True`` only if this
        parameter suffices.

    Returns
    ----------
    object
        First element satisfying this predicate in this iterable if any *or*
        :attr:`betse.util.type.obj.sentinels.SENTINEL` otherwise.

    Raises
    ----------
    BetseIterableException
        If this iterable contains no such item.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.obj.sentinels import SENTINEL

    # Collective efficiency is our middle names.
    return next((item for item in iterable if predicate(item)), SENTINEL)


@type_check
def get_item_first_satisfying(
    iterable: IterableTypes,
    predicate: CallableTypes,
    exception_message: str = None,
) -> object:
    '''
    First item of the passed iterable satisfying the passed **predicate** (i.e.,
    callable accepting one parameter, returning `True` only if this parameter
    suffices) if this iterable contains such an item *or* raise an exception
    otherwise (i.e., if this iterable contains no such item).

    Parameters
    ----------
    iterable : IterableTypes
        Iterable to be inspected.
    predicate : CallableTypes
        Callable accepting one parameter and returning `True` only if this
        parameter suffices.
    exception_message : optional[str]
        Exception message to be raised if no such item is found. Defaults to
        ``None``, in which case a general-purpose message is synthesized.

    Returns
    ----------
    object
        First item satisfying this predicate in this iterable.

    Raises
    ----------
    BetseIterableException
        If this iterable contains no such item.

    See Also
    ----------
    :func:`get_item_first_satisfying_or_sentinel`
        Further details on ordering guarantees.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.obj.sentinels import SENTINEL

    # First item satisfying this predicate in this iterable if any *OR* the
    # sentinel placeholder otherwise.
    first_item = get_item_first_satisfying_or_sentinel(iterable, predicate)

    # If no item satifies this predicate, raise an exception.
    if first_item is SENTINEL:
        # If no exception message is passed, synthesize a default message.
        if exception_message is None:
            exception_message = (
                'Iterable "{}" item satisfying predicate {} not found.'.format(
                    iterable, predicate))

        # Raise this exception.
        raise BetseIterableException(exception_message)

    # Else, return this element.
    return first_item

# ....................{ GETTERS ~ last : instance         }....................
@type_check
def get_item_last_instance_of(
    iterable: IterableTypes, cls: TestableTypes, **kwargs) -> object:
    '''
    Last instance of the passed class or tuple of classes retrieved from the
    passed iterable if this iterable contains such an item *or* raise an
    exception otherwise (i.e., if this iterable contains no such item).

    Parameters
    ----------
    iterable : IterableTypes
        Iterable to be searched.
    cls : TestableTypes
        Class or tuple of classes of the item to search this iterable for.
    kwargs : dict
        Dictionary of all remaining keyword arguments to be passed as is to the
        :func:`get_item_last_satisfying` function.

    Returns
    ----------
    object
        Last instance of this class in this iterable.

    Raises
    ----------
    BetseIterableException
        If this iterable contains no such item.

    See Also
    ----------
    :func:`get_item_last_satisfying`
        Further details on ordering guarantees.
    '''

    return get_item_last_satisfying(
        iterable=iterable,
        predicate=lambda item: isinstance(item, cls),
        **kwargs
    )


@type_check
def get_item_last_instance_of_or_none(
    iterable: IterableTypes, cls: TestableTypes, **kwargs) -> object:
    '''
    Last instance of the passed class or tuple of classes retrieved from the
    passed iterable if this iterable contains such an item *or* ``None``
    otherwise (i.e., if this iterable contains no such element).

    Parameters
    ----------
    iterable : IterableTypes
        Iterable to be searched.
    cls : TestableTypes
        Class or tuple of classes of the item to search this iterable for.
    kwargs : dict
        Dictionary of all remaining keyword arguments to be passed as is to the
        :func:`get_item_last_satisfying` function.

    Returns
    ----------
    object
        Last instance of this class in this iterable if any *or*
        :attr:`betse.util.type.obj.sentinels.SENTINEL` otherwise.

    Raises
    ----------
    BetseIterableException
        If this iterable contains no such item.

    See Also
    ----------
    :func:`get_item_last_satisfying`
        Further details on ordering guarantees.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.obj.sentinels import SENTINEL

    # Last instance of this class in this iterable if any or the sentinel
    # placeholder otherwise.
    item_found = get_item_last_satisfying_or_sentinel(
        iterable=iterable,
        predicate=lambda item: isinstance(item, cls),
        **kwargs
    )

    # Return this instance if *NOT* the sentinal placeholder or None otherwise.
    return item_found if item_found is not SENTINEL else None

# ....................{ GETTERS ~ last : satisfying       }....................
@type_check
def get_item_last_satisfying_or_sentinel(
    iterable: IterableTypes, predicate: CallableTypes) -> object:
    '''
    Last element of the passed iterable satisfying the passed **predicate**
    (i.e., callable accepting one parameter, returning ``True`` only if this
    parameter suffices) if this iterable contains such an element *or* the
    :attr:`betse.util.type.obj.sentinels.SENTINEL` placeholder constant
    otherwise.

    If the passed iterable is a:

    * Sequence, this is guaranteed to be the last such element.
    * Non-sequence (e.g., :class:`set`, :class:`dict`), this should be assumed
      to be a random such element. While most non-sequences guarantee
      predictable order of retrieval assuming no intervening changes, this is a
      fairly unreliable assumption.

    Parameters
    ----------
    iterable : IterableTypes
        Iterable to be inspected.
    predicate : CallableTypes
        Callable accepting one parameter and returning `True` only if this
        parameter suffices.

    Returns
    ----------
    object
        Last element satisfying this predicate in this iterable if any *or*
        :attr:`betse.util.type.obj.sentinels.SENTINEL` otherwise.

    Raises
    ----------
    BetseIterableException
        If this iterable contains no such element.
    '''

    # For simplicity, the existing get_item_first_satisfying_or_sentinel()
    # function is deferred to by returning the first element in the reverse of
    # this iterable satisfying this predicate.
    return get_item_first_satisfying_or_sentinel(
        # For safety, this iterable is reversed via the high-level reverse()
        # function rather than the low-level reversed() builtin; the latter
        # fails to generically support all possible iterable types.
        iterable=reverse(iterable),
        predicate=predicate,
    )



@type_check
def get_item_last_satisfying(
    iterable: IterableTypes,
    predicate: CallableTypes,
    exception_message: str = None,
) -> object:
    '''
    Last element of the passed iterable satisfying the passed **predicate**
    (i.e., callable accepting one parameter, returning `True` only if this
    parameter suffices) if this iterable contains such an element _or_ raise an
    exception otherwise (i.e., if this iterable contains no such element).

    Parameters
    ----------
    iterable : IterableTypes
        Iterable to be inspected.
    predicate : CallableTypes
        Callable accepting one parameter and returning `True` only if this
        parameter suffices.
    exception_message : optional[str]
        Exception message to be raised if no such element is found. Defaults to
        `None`, in which case a suitably general-purpose message is synthesized.

    Returns
    ----------
    object
        Last element satisfying this predicate in this iterable.

    Raises
    ----------
    BetseIterableException
        If this iterable contains no such element.

    See Also
    ----------
    :func:`get_item_first_satisfying_or_sentinel`
        Further details on ordering guarantees.
    '''

    # For simplicity, the existing get_item_first_satisfying() function is
    # deferred to by returning the first element in the reverse of this iterable
    # satisfying this predicate.
    return get_item_first_satisfying(
        # For safety, this iterable is reversed via the high-level reverse()
        # function rather than the low-level reversed() builtin; the latter
        # fails to generically support all possible iterable types.
        iterable=reverse(iterable),
        predicate=predicate,
        exception_message=exception_message,
    )

# ....................{ CONVERTERS                         }....................
@type_check
def to_iterable(iterable: IterableTypes, cls: ClassType) -> IterableTypes:
    '''
    Convert the passed iterable into an iterable of the passed type.

    If this iterable is:

    * Of the same type as the passed type, this iterable is returned as is.
    * A non-Numpy iterable (e.g., :class:`list`) and the passed type is that of:
      * Another non-Numpy iterable (e.g., :class:`tuple`), this iterable is
        converted into an instance of this type. To do so, this type's
        ``__init__`` method is expected to accept this :class:`list` as a single
        positional argument.
      * A Numpy array, this iterable is converted to a Numpy array via the
        :func:`betse.lib.numpy.nparray.from_iterable` function.
    * A Numpy array, this array is converted to the passed type via the
      :func:`betse.lib.numpy.nparray.to_iterable` function.

    Parameters
    ----------
    iterable: IterableTypes
        Source iterable to be converted.
    cls : ClassType
        Type of the target iterable to convert this source iterable into.

    Returns
    ----------
    IterableTypes
        Target iterable converted from this source iterable.
    '''

    # Avoid importing third-party packages at the top level, for safety.
    from betse.lib.numpy import nparray
    from numpy import ndarray

    # Type of the source iterable.
    iterable_src_type = type(iterable)

    # If the source and target iterables are of the same type, return this
    # source iterable as is.
    if iterable_src_type is cls:
        return iterable

    # Else if the source iterable is a Numpy array, defer to logic elsewhere.
    if iterable_src_type is ndarray:
        return nparray.to_iterable(array=iterable, cls=cls)

    # Else if the target iterable is a Numpy array, defer to logic elsewhere.
    if cls is ndarray:
        return nparray.from_iterable(iterable)

    # Else, the source iterable is a non-Numpy iterable. Defer to the
    # constructor of the target iterable for conversion.
    return cls(iterable)

# ....................{ CONSUMERS                          }....................
@type_check
def consume(iterable: IterableTypes, iterations: int) -> object:
    '''
    Consume the passed number of iterations from the passed iterable by
    advancing this iterable forward by this number of iterations.

    For efficiency, this iterable is consumed at C speeds through standard
    objects implemented in low-level C rather than high-level Python.

    Parameters
    ----------
    iterable : IterableTypes
        Iterable to be consumed.
    iterations : int
        Number of iterations to advance this iterable. This number should be
        strictly positive (i.e., `iterations >= 1`).

    Returns
    ----------
    object
        Object yielded by the last iterable iteration (i.e., the last `next()`
        method called on this iterable) if any _or_ `None` if this iterable was
        already exhausted (i.e., empty) when passed.

    See Also
    ----------
    http://docs.python.org/3/library/itertools.html
        Official documentation strongly inspiring this function.
    '''
    assert types.is_int_positive(iterations), (
        types.assert_not_int_positive(iterations))

    # Advance to the empty slice starting (and hence immediately ending) at the
    # passed iterable position and return the value of this iterable at this
    # position if any.
    return next(slice(iterable, iterations, iterations), None)


@type_check
def exhaust(iterable: IterableTypes) -> object:
    '''
    Exhaust the passed iterable by advancing this iterable directly past its
    last iteration.

    For efficiency, this iterable is consumed at C speeds through standard
    objects implemented in low-level C rather than high-level Python.

    Caveats
    ----------
    **This function should only be called for finite iterables.** If the passed
    iterable:

    * Explicitly halts with a :class:`StopIteration` exception and hence is
      finite, this function exhausts this iterable as expected.
    * Does *not* explicitly halt with a `StopIteration` exception and hence is
      infinite, this function reduces to an **infinite loop.**

    Parameters
    ----------
    iterable : IterableTypes
        Iterable to be exhausted.

    Returns
    ----------
    object
        Object yielded by the last iterable iteration (i.e., the last `next()`
        method called on this iterable) if any _or_ `None` if this iterable was
        already exhausted (i.e., empty) when passed.

    See Also
    ----------
    http://docs.python.org/3/library/itertools.html
        Official documentation strongly inspiring this function.
    '''

    # For efficiency, feed this iterable into a zero-length deque retaining only
    # the last iterated value if any.
    iterable_deque = deque(iterable, maxlen=1)

    # If this iterable was *NOT* already exhausted, return its last value.
    if iterable_deque:
        return iterable_deque[0]
    # Else, this iterable was already exhausted. Return nothing.
    else:
        return None

# ....................{ INVERTERS                          }....................
@type_check
def invert_iterable_unique(iterable: IterableTypes) -> MappingType:
    '''
    Dictionary inverted from the passed iterable if **internally unique** (i.e.,
    containing no duplicate items) *or* raise an exception otherwise.

    Specifically:

    * If this iterable is a dictionary, the
      :func:`betse.util.type.mapping.mappings.invert_dict_unique` function is
      silently deferred to. The type of the returned dictionary is guaranteed to
      be the same as the type of the passed dictionary.
    * Else, the returned dictionary maps from each item of this iterable to the
      0-based index of that item in this iterable. The type of this dictionary
      is guaranteed to *not* be the same as the type of this iterable.

    Parameters
    ----------
    iterable : IterableTypes
        Internally unique iterable to be inverted.

    Returns
    ----------
    MappingType
        Dictionary inverted from this iterable as detailed above.

    Raises
    ----------
    BetseIterableException
        If at least one item of this iterable is a duplicate.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.mapping import mappings

    # If this iterable is a mapping...
    if mappings.is_mapping(iterable):
        # Invert this iterable specifically as a mapping.
        #
        # While mappings are technically iterables and hence invertable via the
        # generic approach applied below, doing so incorrectly returns a
        # dictionary mapping from the keys of the passed dictionary to arbitrary
        # 0-based integers. As the adjective "arbitrary" implies, this renders
        # the resulting dictionary effectively useless for most purposes.
        return mappings.invert_map_unique(iterable)
    # Else, this iterable is a non-mapping. In this case, a generic approach
    # suffices... usually.
    else:
        # If one or more items of this iterable are duplicates, raise an
        # exception.
        die_unless_items_unique(iterable)

        # One-liners for Great Glory.
        return {item: item_index for item_index, item in enumerate(iterable)}

# ....................{ ITERATORS                          }....................
@type_check
def iter_items(*iterables: IterableTypes) -> GeneratorType:
    '''
    Generator yielding each item of each of the passed iterables (in both the
    internal order of each iterable and the passed order of iterables),
    effectively "chaining" these iterables together.

    This function is preferable for pure iteration over multiple iterables, in
    which case a composite iterable of the same type is *not* required.

    Parameters
    ----------
    iterables : tuple[IterableTypes]
        Tuple of all iterables whose items are to be iterated over.

    Yields
    ----------
    object
        Current item of the current iterable.

    See Also
    ----------
    :func:`join`
        Less efficient alternative producing a composite iterable.
    https://stackoverflow.com/a/14342360/2809027
        Stackoverflow answer strongly inspiring this implementation.

    Examples
    ----------
        >>> from betse.util.type import iterables
        >>> in_xanadu = ('did', 'Kubla', 'Khan')
        >>> a_stately = ('pleasure-dome', 'decree:')
        >>> for word in iterables.iter_items(in_xanadu, a_stately):
        ...     print(word, end=' ')
        did Kubla Khan pleasure-dome decree:
    '''

    # Pure Python 3.x. Join the party.
    #
    # Note that the itertools.chain.from_iterable() could also be called here,
    # but that doing so is less Pythonic than the current approach.
    for iterable in iterables:
        yield from iterable

# ....................{ JOINERS                            }....................
@type_check
def join(*iterables: IterableTypes) -> IterableTypes:
    '''
    Join each item of each of the passed iterables (in both the internal order
    of each iterable and the passed order of iterables) into a new iterable of
    the same type as the first such iterable.

    Parameters
    ----------
    iterables : tuple[IterableTypes]
        Tuple of all iterables whose items are to be joined together.

    Returns
    ----------
    IterableTypes
        Iterable of the same type as the first passed iterable containing each
        item of each of the passed iterables.

    Raises
    ----------
    BetseIterableException
        If the type of the first passed iterable is :class:`str`, in which case
        the builtin :func:`str.join` function should be called for both
        efficiency and sanity instead.

    See Also
    ----------
    :func:`iter_items`
        More efficient alternative intended for pure iteration, in which case a
        composite iterable of the same type is *not* required.
    https://stackoverflow.com/a/14342360/2809027
        Stackoverflow answer strongly inspiring this implementation.

    Examples
    ----------
        >>> from betse.util.type import iterables
        >>> where_alph = ('the', 'sacred', 'river,', 'ran')
        >>> through_caverns = ('measureless', 'to', 'man')
        >>> for word in iterables.iter_items(where_alph, through_caverns):
        ...     print(word, end=' ')
        the sacred river, ran measureless to man
    '''

    # Iterator over each of the passed iterables.
    iterator = iter(iterables)

    # First passed iterable.
    iterable_first = next(iterator)

    # Type of this iterable.
    iterable_first_type = type(iterable_first)

    # If this iterable is a string, raise an exception.
    if iterable_first_type is str:
        raise BetseIterableException(
            'String "{}" not joinable by iterables.join(). '
            'Consider calling str.join() instead.'.format(iterable_first))

    # Return a new iterable of this type over each item of each passed iterable.
    return iterable_first_type(iter_items(*iterables))

# ....................{ REVERSERS                          }....................
@type_check
def reverse(iterable: IterableTypes) -> IterableTypes:
    '''
    Reverse the passed iterable into a new iterable of differing type containing
    all elements of the passed iterable in reverse order.

    Parameters
    ----------
    iterable : IterableTypes
        Iterable to be returned reversed. For generality, this iterable is _not_
        modified by this function.

    Returns
    ----------
    IterableTypes
        Iterable reversed from the passed iterable. For efficiency, this
        iterable is only a shallow rather than deep copy of the passed iterable.
    '''

    # If this iterable is *NOT* reversible as is, convert this iterable into the
    # most space- and time-efficient iterable containing the same elements that
    # *IS* reversible -- in this case, a tuple.
    if not is_reversible(iterable):
        iterable = tuple(iterable)

    # Return the result of the efficient reversed() builtin on this iterable,
    # now guaranteed to be reversible as is.
    return reversed(iterable)

# ....................{ SORTERS ~ ascending                }....................
@type_check
def sort_ascending(iterable: IterableTypes) -> IterableTypes:
    '''
    Iterable sorted from the passed iterable in ascending order.

    Each element of this iterable is compared to each other element of this
    iterable via the `<` operator, implicitly calling the `__le__()` special
    method of these elements. Each element is ideally but _not_ necessarily of
    the same type. If each element is:

    * A string, these strings are sorted in **ascending lexicographic order**
      (i.e., traditional order of dead-tree dictionaries and encyclopedias).
    * A number (i.e., integer or a float), these numbers are sorted in
      **ascending numeric order.**

    Parameters
    ----------
    iterable : IterableTypes
        Unsorted iterable to be returned sorted. For generality, this iterable
        is _not_ modified by this function.

    Returns
    ----------
    IterableTypes
        Iterable sorted from and of the same type as the passed iterable. For
        efficiency, this iterable is only a shallow rather than deep copy of the
        passed iterable. Note lastly that the class of the passed iterable
        _must_ define an `__init__()` method accepting a list.
    '''

    # Type of the passed iterable.
    iterable_type = type(iterable)

    # Return an iterable of the same type, converted from the sorted list
    # returned by the sorted() builtin.
    return iterable_type(sorted(iterable))


@type_check
def sort_by_index_ascending(
    iterable: IterableTypes, subiterable_index: object) -> IterableTypes:
    '''
    Iterable of subiterables sorted from the passed iterable of subiterables in
    ascending order of the value of each element at the passed key or index of
    each subiterable of this iterable.

    Each element at the passed key or index of each subiterable of this iterable
    is compared to each other element at each other key or index of each other
    subiterable of this iterable via the `<` operator, implicitly calling the
    `__le__()` special method of these elements. Each element is ideally but
    _not_ necessarily of the same type. If each element is:

    * A string, these strings are sorted in **ascending lexicographic order**
      (i.e., traditional order of dead-tree dictionaries and encyclopedias).
    * A number (i.e., integer or a float), these numbers are sorted in
      **ascending numeric order.**

    Parameters
    ----------
    iterable : IterableTypes
        Unsorted iterable of subiterables to be returned sorted. For
        generality, neither this iterable nor these subiterables are modified
        by this function.
    subiterable_index : object
        Object with which to index each subiterable of this iterable. The type
        of this object _must_ be a type accepted by the `__getitem__()` special
        method of each subiterable. Specifically, if each subiterable is a:
        * **Mapping** (e.g., :class:`dict`), this object _must_ be hashable.
        * **Sequence** (e.g., :class:`list`, :class:`tuple`), this object
          _must_ be either:
          * An integer.
          * A :func:`slice` object.

    Returns
    ----------
    IterableTypes
        Iterable of subiterables sorted from and of the same type as the passed
        iterable of subiterables. For efficiency, this iterable is only a
        shallow rather than deep copy of the passed iterable. Note lastly that
        the class of the passed iterable _must_ define an `__init__()` method
        accepting a list.
    '''

    # Type of the passed iterable.
    iterable_type = type(iterable)

    # Return an iterable of the same type, converted from the sorted list
    # returned by the sorted() builtin.
    #
    # For efficiency, elements of each subiterable of this iterable are
    # retrieved via the function internally created and called by calling an
    # instance of the standard "itemgetter" class. While technically an
    # instance of this class, this class is intended to be treated as a simple
    # function passed the desired index returning a function passed a
    # subiterable returning the element at that index: e.g.,
    #
    #     >>> penguins = [('adelie', 0xFEEDFACE), ('gentoo', 0xDEADBEEF)]
    #     >>> penguin_id = itemgetter(2)
    #     >>> penguin_id(penguins[0]) == 0xFEEDFACE
    #     True
    #
    # Despite the internal complexity and object overhead imposed by the
    # "itemgetter" class, this approach has been definitively profiled to be
    # 126% faster on average than the traditional
    # "lambda subiterable: subiterable[subiterable_index]" approach --
    # presumably due to hidden overhead imposed by capturing all local
    # variables into the closure context encapsulated by the lambda. See also
    # the following Stackoverflow answer exhibiting this profiling:
    #
    #     https://stackoverflow.com/a/17243726/2809027
    return iterable_type(sorted(iterable, key=itemgetter(subiterable_index)))

# ....................{ SORTERS ~ descending               }....................
@type_check
def sort_descending(iterable: IterableTypes) -> IterableTypes:
    '''
    Iterable sorted from the passed iterable in descending order.

    Each element of this iterable is compared to each other element of this
    iterable via the `>` operator, implicitly calling the `__ge__()` special
    method of these elements. Each element is ideally but _not_ necessarily of
    the same type. If each element is:

    * A string, these strings are sorted in **descending lexicographic order**
      (i.e., reverse order of dead-tree dictionaries and encyclopedias).
    * A number (i.e., integer or a float), these numbers are sorted in
      **descending numeric order.**

    See Also
    ----------
    :func:`sort_ascending`
        Further details.
    '''

    # Type of the passed iterable.
    iterable_type = type(iterable)

    # Return an iterable of the same type, converted from the sorted list
    # returned by the sorted() builtin.
    return iterable_type(sorted(iterable, reverse=True))


@type_check
def sort_by_index_descending(
    iterable: IterableTypes, subiterable_index: object) -> IterableTypes:
    '''
    Iterable of subiterables sorted from the passed iterable of subiterables in
    descending order of the value of each element at the passed key or index of
    each subiterable of this iterable.

    Each element at the passed key or index of each subiterable of this
    iterable is compared to each other element at each other key or index of
    each other subiterable of this iterable via the `>` operator, implicitly
    calling the `__ge__()` special method of these elements. Each element is
    ideally but _not_ necessarily of the same type. If each element is:

    * A string, these strings are sorted in **descending lexicographic order**
      (i.e., reverse order of dead-tree dictionaries and encyclopedias).
    * A number (i.e., integer or a float), these numbers are sorted in
      **descending numeric order.**

    See Also
    ----------
    :func:`sort_by_index_ascending`
        Further details.
    '''

    # Type of the passed iterable.
    iterable_type = type(iterable)

    # Return an iterable of the same type, converted from the sorted list
    # returned by the sorted() builtin. See sort_by_index_descending() for
    # commentary, particularly on the "key" parameter.
    return iterable_type(sorted(
        iterable, key=itemgetter(subiterable_index), reverse=True))

# ....................{ ZIPPERS                            }....................
#FIXME: Unit test us up.
@type_check
def zip_isometric(*iterables: IterableTypes) -> GeneratorType:
    '''
    Generator zipping all passed iterables required to be of the same length.

    This generator iteratively yields an `n`-tuple, where:

    * `n` is the length of each passed iterable.
    * The `i`-th element of this tuple is in the `i`-th passed iterable.

    Parameters
    ----------
    iterables : IterableTypes
        Tuple of iterables of the same length to be zipped.

    Returns
    ----------
    GeneratorType
        Generator zipping these iterables.

    Raises
    ----------
    BetseIterableException
        If any passed iterable differs in length from any other passed iterable.

    See Also
    ----------
    https://stackoverflow.com/a/32954700/2809027
        Stackoverflow answer strongly inspiring this implementation.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.obj.sentinels import SENTINEL

    # Iteratively zip and yield each n-tuple from the passed n iterables. To
    # efficiently detect iterables of insufficient length, the C-based
    # zip_longest() function is called to fill all iterables of insufficient
    # length with sentinel objects to the expected length.
    #
    # After zipping but before yielding each n-tuple, this n-tuple is then
    # manually searched for sentinel objects. Since these objects may reside at
    # any index of this n-tuple, the entire n-tuple *MUST* be searched in an
    # O(n) manner. While unfortunate, this approach remains substantially more
    # efficient than all alternatives -- largely due to the efficacy of the
    # zip_longest() function. Surprisingly, this "for" loop-based approach has
    # been timed to be faster than the following generator expression:
    #
    #     return (
    #         ntuple if SENTINEL not in ntuple else (
    #             _zip_isometric_error(iterables, ntuple))
    #         for ntuple in zip_longest(
    #             *iterables, fillvalue=SENTINEL)
    #     )
    for ntuple in itertools.zip_longest(
        *iterables, fillvalue=SENTINEL):
        # If this n-tuple contains a sentinel, at least one passed iterable is
        # of insufficient length. Raise a human-readable exception indicating
        # the index and contents of this iterable.
        if SENTINEL in ntuple:
            raise _zip_isometric_error(iterables, ntuple)

        # Else, this n-tuple is valid. Yield it up!
        yield ntuple


def _zip_isometric_error(iterables: tuple, ntuple: tuple) -> None:
    '''
    Raise an exception indicating that the iterable of the passed tuple of
    iterables identified by the passed `n`-tuple is of smaller length than other
    iterables in this tuple of iterables.

    This private function is _only_ intended to be called by the
    :func:`zip_isometric` function.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.obj.sentinels import SENTINEL

    # Index of the erroneously short iterable in this tuple of iterables,
    # identical to the index of the first sentinel in the passed zipped tuple.
    iterable_short_index = None
    for iterable_short_index, item in enumerate(ntuple):
        if item is SENTINEL:
            break

    # This erroneously short iterable.
    iterable_short = iterables[iterable_short_index]

    # If the first iterable is of predefined length (e.g., is *NOT* a generator
    # of dynamic length), end this exception message with this length.
    if isinstance(iterables[0], SizedType):
        exception_suffix = 'length {} of prior iterables'.format(
            len(iterables[0]))
    # Else, end this exception message as is.
    else:
        exception_suffix = 'length of prior iterables'

    # If this erroneously short iterable is of predefined length (e.g., is *NOT*
    # a generator of dynamic length), begin this exception message with this
    # length and end this message with this iterable's contents.
    if isinstance(iterable_short, SizedType):
        exception_prefix = "Length {} of iterable {}".format(
            len(iterable_short), iterable_short_index)
        exception_suffix += ': {!r}'.format(iterable_short)
    # Else, begin this exception message as is.
    else:
        exception_prefix = "Length of iterable {}".format(iterable_short_index)

    # Raise this exception.
    raise BetseIterableException('{} differs from {}'.format(
        exception_prefix, exception_suffix))
