#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level **non-string iterable getters** (i.e., utility functions getting one
or more items contained within non-string objects implementing the abstract
:class:`collections.abc.Iterable` base class) facilities.
'''

# ....................{ IMPORTS                           }....................
from betse.exceptions import BetseIterableException
from betse.util.io.log import logs
from betse.util.type.types import (
    type_check, CallableTypes, IterableTypes, TestableTypes)

# ....................{ GETTERS                           }....................
@type_check
def get_items_duplicate(iterable: IterableTypes) -> set:
    '''
    Unordered set of all duplicate items in the passed iterable.

    Parameters
    ----------
    iterable : IterableTypes
        Iterable to be inspected.

    Returns
    ----------
    set
        Unordered set of all duplicate items in this iterable.

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

    # Return the set of all duplicate items.
    return items_duplicate

# ....................{ GETTERS ~ first                   }....................
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

# ....................{ GETTERS ~ first : instance        }....................
@type_check
def get_item_first_instance_of(
    iterable: IterableTypes, cls: TestableTypes, **kwargs) -> object:
    '''
    First instance of the passed class or tuple of classes retrieved from the
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
    First item of this iterable that is *not* an instance of the passed class
    or tuple of classes if this iterable contains such an item *or* raise an
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

# ....................{ GETTERS ~ first : satisfying      }....................
@type_check
def get_item_first_satisfying_or_sentinel(
    iterable: IterableTypes, predicate: CallableTypes) -> object:
    '''
    First item of the passed iterable satisfying the passed **predicate**
    (i.e., callable accepting one parameter, returning ``True`` only if this
    parameter suffices) if this iterable contains such an item *or* the
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
    First item of the passed iterable satisfying the passed **predicate**
    (i.e., callable accepting one parameter, returning ``True`` only if this
    parameter suffices) if this iterable contains such an item *or* raise an
    exception otherwise (i.e., if this iterable contains no such item).

    Parameters
    ----------
    iterable : IterableTypes
        Iterable to be inspected.
    predicate : CallableTypes
        Callable accepting one parameter and returning ``True`` only if this
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
        Callable accepting one parameter and returning ``True`` only if this
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

    # Avoid circular import dependencies.
    from betse.util.type.iterable import iterables

    # For simplicity, the existing get_item_first_satisfying_or_sentinel()
    # function is deferred to by returning the first element in the reverse of
    # this iterable satisfying this predicate.
    return get_item_first_satisfying_or_sentinel(
        # For safety, this iterable is reversed via the high-level reverse()
        # function rather than the low-level reversed() builtin; the latter
        # fails to generically support all possible iterable types.
        iterable=iterables.reverse(iterable),
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
    (i.e., callable accepting one parameter, returning ``True`` only if this
    parameter suffices) if this iterable contains such an element *or* raise an
    exception otherwise (i.e., if this iterable contains no such element).

    Parameters
    ----------
    iterable : IterableTypes
        Iterable to be inspected.
    predicate : CallableTypes
        Callable accepting one parameter and returning ``True`` only if this
        parameter suffices.
    exception_message : optional[str]
        Exception message to be raised if no such element is found. Defaults to
        ``None``, in which case a suitably general-purpose message is
        synthesized.

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

    # Avoid circular import dependencies.
    from betse.util.type.iterable import iterables

    # For simplicity, the existing get_item_first_satisfying() function is
    # deferred to by returning the first element in the reverse of this
    # iterable satisfying this predicate.
    return get_item_first_satisfying(
        # For safety, this iterable is reversed via the high-level reverse()
        # function rather than the low-level reversed() builtin; the latter
        # fails to generically support all possible iterable types.
        iterable=iterables.reverse(iterable),
        predicate=predicate,
        exception_message=exception_message,
    )

# ....................{ GETTERS ~ str                     }....................
@type_check
def get_item_var_uniquified_str(
    iterable: IterableTypes,
    item_var_name: str,
    item_var_format: str,
) -> str:
    '''
    Create and return a new machine-readable string guaranteed to be unique
    across each instance variable with the passed name of each item of the
    passed iterable, typically to enforce uniqueness of an instance variable
    employed by the caller as a **primary key** (i.e., attribute uniquely
    identifying each item of this iterable).

    This function requires each item of this iterable to declare an attribute
    with the passed name whose value is an arbitrary string. This function then
    synthesizes a string suitable for use by callers as the value of that
    attribute for a currently non-existing item of this iterable (to be
    presumably created by the caller after calling this function).

    Parameters
    ----------
    iterable : IterableTypes
        Iterable to be inspected.
    item_var_name : str
        Name of the instance variable declared by *all* items of this list,
        whose string values are to be uniquified.
    item_var_format : str
        Format specifier containing a ``{}`` substring (e.g., ``Item ({}).``),
        iteratively interpolated by this function with an arbitrary integer to
        produce the returned string.

    Returns
    ----------
    str
        New machine-readable string guaranteed to both match this format *and*
        be unique across all instance variables with this name of all items of
        this iterable.

    Raises
    ----------
    BetseObjectException
        If, for some item of this iterable, either:

        * This item contains no instance variable with this name.
        * This item contains an instance variable with this name whose value is
          *not* a string.
    BetseStrException
        If the passed format specifier contains no ``{}`` substring.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.obj import objects
    from betse.util.type.text.string import strs

    # Log this formatting.
    logs.log_debug(
        'Uniquifying iterable item variable "%s" with template "%s"...',
        item_var_name, item_var_format)

    #FIXME: Ideally, we would also raise exceptions if this format
    #specifier contains two or more ``{}`` substrings. Sadly, there appears
    #no trivial and/or efficient means of doing so. While we suppose we
    #*COULD* define a new
    #betse.util.type.text.string.strs.get_substrs_count() function doing
    #so, the lack of any benefit to doing so hardly seems worth it. Note
    #also that genuinely testing this is effectively infeasible due to
    #false positives (e.g., the string "{}{}" should raise exceptions but
    #the string "{}{{}}" should *NOT*).

    # If this specifier contains no "{}" substring, raise an exception.
    #
    # Note that this simplistic logic fails to account for "{{" and "}}"
    # escaping and hence *COULD* fail to raise exceptions when passed
    # worst-case format specifiers, but that we mostly do not care.
    strs.die_unless_substr(text=item_var_format, substr='{}')

    # Arbitrary unique identifier with which to uniquify (i.e., guarantee
    # the uniqueness of) an arbitrary string attribute of a new item in
    # this list, defaulting to the number of existing items in this list.
    item_id = len(iterable)

    # Current string value of the instance variable with this name of the
    # current item of this iterable.
    item_var = None

    # Uniquified string to be created and returned by this function.
    item_var_new = None

    # While this string is *NOT* unique across this iterable, iteratively
    # (re)search this iterable until obtaining a unique string. Assuming this
    # iterable does *not* define special methods in non-sane ways, iteration
    # is guaranteed to terminate successfully.
    while True:
        item_var_new = item_var_format.format(item_id)

        # For each existing item of this iterable...
        for item in iterable:
            # Current string value of the instance variable with this name of
            # the current item of this iterable if this item defines such a
            # variable *OR* raise an exception otherwise.
            item_var = objects.get_attr(
                obj=item, attr_name=item_var_name, attr_type=str)

            # If this value collides with the uniquified string to be created
            # and returned by this function, this string is non-unique. In this
            # case, increment this identifier, format a new attribute with this
            # identifier, and continue searching.
            if item_var == item_var_new:
                item_id += 1
                break
        # If the above iteration terminated successfully, this string is
        # unique. In this case, halt searching and return this string.
        else:
            break
        # Else, the above iteration terminated prematurely, implying this
        # string to be non-unique. In this case, continue searching.

    # Return this string.
    return item_var_new
