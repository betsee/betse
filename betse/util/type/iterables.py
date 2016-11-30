#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry
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
    GeneratorType,
    IterableTypes,
    SizedType,
    TestableTypes,
)
from collections import deque
from operator import itemgetter

# ....................{ CLASSES                            }....................
class Sentinel(object):
    '''
    Class encapsulating sentinel objects of arbitrary (albeit unique) value.

    Instances of this class are intended to be used as placeholder objects in
    iterables, typically to identify erroneous or edge-case algorithm input.
    '''

    def __repr__(self) -> str:
        '''
        Human- and machine-readable representation of this sentinel.

        This method has been overridden purely to improve the debuggability of
        algorithms requiring instances of this class.
        '''

        return 'Sentinel()'

# ....................{ CONSTANTS                          }....................
SENTINEL = Sentinel()
'''
Sentinel object of arbitrary value.

This object is internally leveraged by various utility functions (e.g.,
:func:`zip_isometric`) to identify erroneous and edge-case input (e.g.,
iterables of insufficient length).
'''

# ....................{ TESTERS                            }....................
@type_check
def is_reversible(iterable: IterableTypes) -> bool:
    '''
    `True` only if the passed iterable is **reversible** (i.e., successfully
    passable to the :func:`reversed` builtin).

    This function returns `True` if this iterable either:

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
    from betse.util.type import objects

    # Return True only if this iterable either...
    return (
        # Defines the __reversed__() special method.
        objects.is_method(iterable, '__reversed__') or (
            # Defines the __len__() and __getitem__() special methods.
            objects.is_method(iterable, '__len__') and
            objects.is_method(iterable, '__getitem__')
        )
    )

# ....................{ TESTERS ~ item                     }....................
@type_check
def is_items_any_satisfying(
    iterable: IterableTypes, predicate: CallableTypes) -> bool:
    '''
    `True` only if some element of the passed iterable satisfies the passed
    **predicate** (i.e., callable accepting one parameter, returning `True` only
    if this parameter suffices).

    Parameters
    ----------
    iterable : IterableTypes
        Iterable to be inspected.
    predicate : CallableTypes
        Callable accepting one parameter and returning `True` only if this
        parameter suffices.

    Returns
    ----------
    bool
        `True` only if some element of this iterable satisfies this predicate.
    '''

    # First item satisfying this predicate in this iterable if any *OR* the
    # sentinel placeholder otherwise.
    first_item = get_item_first_satisfying_or_sentinel(iterable, predicate)

    # Return True only if this item is *NOT* the sentinel, in which case some
    # item satisfies this predicate.
    return first_item is not SENTINEL


@type_check
def is_items_any_instance_of(
    iterable: IterableTypes, cls: TestableTypes) -> bool:
    '''
    `True` only if some element of the passed iterable is an instance of the
    passed class.

    Parameters
    ----------
    iterable : IterableTypes
        Iterable to be inspected.
    cls : TestableTypes
        Type of the element to be tested for.

    Returns
    ----------
    bool
        `True` only if some element of this iterable is an instance of this
        class.
    '''

    return is_items_any_satisfying(
        iterable=iterable, predicate=lambda item: isinstance(item, cls))

# ....................{ GETTERS ~ first                    }....................
@type_check
def get_item_first(iterable: IterableTypes) -> object:
    '''
    First element non-destructively retrieved from the passed iterable if this
    iterable is non-empty _or_ raise an exception otherwise (i.e., if this
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
        Cecil's Stackoverflow answer strongly inspiring this function, complete
        with detailed timings of all alternative solutions.
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
    for first_item in iterable:
        break

    # Return the first element iterated above.
    return first_item


@type_check
def get_item_first_satisfying_or_sentinel(
    iterable: IterableTypes, predicate: CallableTypes) -> object:
    '''
    First element of the passed iterable satisfying the passed **predicate**
    (i.e., callable accepting one parameter, returning `True` only if this
    parameter suffices) if this iterable contains such an element _or_ the
    :data:`SENTINEL` placeholder constant otherwise.

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
        Callable accepting one parameter and returning `True` only if this
        parameter suffices.

    Returns
    ----------
    object
        First element satisfying this predicate in this iterable if any _or_
        :data:`SENTINEL` otherwise.

    Raises
    ----------
    BetseIterableException
        If this iterable contains no such element.
    '''

    # Collective efficiency is our middle names.
    return next((item for item in iterable if predicate(item)), SENTINEL)


@type_check
def get_item_first_satisfying(
    iterable: IterableTypes,
    predicate: CallableTypes,
    exception_message: str = None,
) -> object:
    '''
    First element of the passed iterable satisfying the passed **predicate**
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
        First element satisfying this predicate in this iterable.

    Raises
    ----------
    BetseIterableException
        If this iterable contains no such element.

    See Also
    ----------
    :func:`get_item_first_satisfying_or_sentinel`
        Further details on ordering guarantees.
    '''

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


@type_check
def get_item_first_instance_of(
    iterable: IterableTypes, cls: TestableTypes, **kwargs) -> object:
    '''
    First instance of the passed class retrieved from the passed iterable if
    this iterable such an element _or_ raise an exception otherwise (i.e., if
    this iterable contains no such element).

    Parameters
    ----------
    iterable : IterableTypes
        Iterable to be inspected.
    cls : TestableTypes
        Type of the element to be retrieved.
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
        If this iterable contains no such element.

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

# ....................{ GETTERS ~ last                    }....................
@type_check
def get_item_last_satisfying_or_sentinel(
    iterable: IterableTypes, predicate: CallableTypes) -> object:
    '''
    Last element of the passed iterable satisfying the passed **predicate**
    (i.e., callable accepting one parameter, returning `True` only if this
    parameter suffices) if this iterable contains such an element _or_ the
    :data:`SENTINEL` placeholder constant otherwise.

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
        Last element satisfying this predicate in this iterable if any _or_
        :data:`SENTINEL` otherwise.

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


@type_check
def get_item_last_instance_of_or_sentinel(
    iterable: IterableTypes, cls: TestableTypes, **kwargs) -> object:
    '''
    Last instance of the passed class retrieved from the passed iterable if
    this iterable such an element _or_  the :data:`SENTINEL` placeholder
    constant otherwise (i.e., if this iterable contains no such element).

    Parameters
    ----------
    iterable : IterableTypes
        Iterable to be inspected.
    cls : TestableTypes
        Type of the element to be retrieved.
    kwargs : dict
        Dictionary of all remaining keyword arguments to be passed as is to the
        :func:`get_item_last_satisfying` function.

    Returns
    ----------
    object
        Last instance of this class in this iterable if any _or_
        :data:`SENTINEL` otherwise.

    Raises
    ----------
    BetseIterableException
        If this iterable contains no such element.

    See Also
    ----------
    :func:`get_item_last_satisfying`
        Further details on ordering guarantees.
    '''

    return get_item_last_satisfying_or_sentinel(
        iterable=iterable,
        predicate=lambda item: isinstance(item, cls),
        **kwargs
    )


@type_check
def get_item_last_instance_of(
    iterable: IterableTypes, cls: TestableTypes, **kwargs) -> object:
    '''
    Last instance of the passed class retrieved from the passed iterable if
    this iterable such an element _or_ raise an exception otherwise (i.e., if
    this iterable contains no such element).

    Parameters
    ----------
    iterable : IterableTypes
        Iterable to be inspected.
    cls : TestableTypes
        Type of the element to be retrieved.
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
        If this iterable contains no such element.

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

# ....................{ CONSUMERS                          }....................
@type_check
def consume(iterable: IterableTypes, iterations: int) -> object:
    '''
    Consume the passed number of iterations from the passed iterable by
    advancing this iterable forward by this number of iterations.

    For efficiency, this iterable is consumed at C speeds via standard objects
    implemented in low-level C rather than high-level Python.

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

    For efficiency, this iterable is exhausted at C speeds via standard objects
    implemented in low-level C rather than high-level Python.

    Caveats
    ----------
    **This function should _only_ be called for finite iterables.** If the
    passed iterable:

    * Explicitly halts with a `StopIteration` exception and hence is finite,
      this function exhausts this iterable as expected.
    * Does _not_ explicitly halt with a `StopIteration` exception and hence is
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

# ....................{ OPERATORS                          }....................
@type_check
def sum_by_index(iterable: IterableTypes, subiterable_index: object) -> object:
    '''
    Summation of each element at the passed key or index of each subiterable of
    the passed iterable.

    Each element at the passed key or index of each subiterable of this
    iterable is totalized via the `+` operator, implicitly calling the
    `__sum__()` special method of these elements. Each element is ideally but
    _not_ necessarily of the same type. If each element is:

    * A string, these strings are concatenated into a single string.
    * A number (i.e., integer or float), these numbers are summed to a single
      number of the widest type of these numbers. Specifically:
      * If any such number is a float, the returned number is also a float.
      * Else, the returned number is an integer.

    Parameters
    ----------
    iterable : IterableTypes
        Iterable of subiterables to be summed.
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
    object
        Object produced by summing each element at this key or index of each
        subiterable of this iterable.
    '''

    # Efficiency and simplicity combine here to form MegaFastSimple.
    return sum(subiterable[subiterable_index] for subiterable in iterable)

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
    :exc:`BetseIterableException`
        If any passed iterable differs in length from any other passed iterable.

    See Also
    ----------
    https://stackoverflow.com/a/32954700/2809027
        Stackoverflow answer strongly inspiring this implementation.
    '''

    # For efficiency, declare this frequently accessed variable to be global,
    # preventing Python from continually attempting to access this variable as a
    # local in the loop below. (This has been timed.)
    global SENTINEL

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

    # Index of the erroneously short iterable in this tuple of iterables,
    # identical to the index of the first sentinel in the passed zipped tuple.
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
