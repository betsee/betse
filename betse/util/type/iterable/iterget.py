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
from betse.exceptions import BetseIterableException, BetseParamException
from betse.util.io.log import logs
from betse.util.type.types import (
    type_check,
    CallableTypes,
    IterableTypes,
    StrOrNoneTypes,
    TestableTypes,
)

# ....................{ GETTERS                           }....................
@type_check
def get_items_duplicate(iterable: IterableTypes) -> set:
    '''
    Unordered set of all duplicate items of the passed iterable.

    Parameters
    ----------
    iterable : IterableTypes
        Iterable to be inspected.

    Returns
    ----------
    set
        Unordered set of all duplicate items of this iterable.

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
def get_item_first(*iterables: IterableTypes) -> object:
    '''
    First item non-destructively retrieved from the first passed non-empty
    iterable if at least one such iterable is non-empty *or* raise an exception
    otherwise (i.e., if all passed iterables are empty).

    Parameters
    ----------
    iterables : tuple[IterableTypes]
        Tuple of all iterables to be inspected.

    Returns
    ----------
    object
        First item of the first passed non-empty iterable.

    Raises
    ----------
    BetseIterableException
        If all passed iterables are empty.

    See Also
    ----------
    :func:`get_item_first_or_sentinel`
        Further details.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.obj.sentinels import SENTINEL

    # First item of these iterables if any *OR* the sentinel placeholder,
    item_first = get_item_first_or_sentinel(*iterables)

    # If no such item exists, raise an exception.
    if item_first is SENTINEL:
        raise BetseIterableException('Iterables empty.')
    # Else, this item exists.

    # Return this item.
    return item_first


@type_check
def get_item_first_or_sentinel(*iterables: IterableTypes) -> object:
    '''
    First item non-destructively retrieved from the first passed non-empty
    iterable if at least one such iterable is non-empty *or* the sentinel
    placeholder otherwise (i.e., if all passed iterables are empty).

    Specifically, if the first passed non-empty iterable is a:

    * Sequence (e.g., :class:`list`, :class:`tuple`), this function returns the
      first item of this sequence.
    * Non-sequence (e.g., :class:`set`, :class:`dict`), this function returns a
      pseudo-random item of this non-sequence. While most non-sequences
      guarantee predictable order of retrieval assuming no intervening changes,
      this is a fairly unreliable assumption.

    Parameters
    ----------
    iterables : tuple[IterableTypes]
        Tuple of all iterables to be inspected.

    Returns
    ----------
    object
        Either:

        * If at least one such iterable is non-empty, the first item of this
          iterable.
        * Else, the sentinel singleton.

    See Also
    ----------
    https://stackoverflow.com/a/40054478/2809027
        Stackoverflow answer strongly inspiring this implementation, complete
        with detailed timings of all alternative solutions.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.obj.sentinels import SENTINEL
    from betse.util.type.iterable import iterables as itermod

    # First item of these iterables if any *OR* the sentinel placeholder,
    # defaulting to the latter for implementation simplicity.
    item_first = SENTINEL

    # If at least one such iterable is non-empty, iterate this iterable once
    # and assign the first item of this iterable to this local when doing so.
    #
    # Note that this implementation, despite being verbose and unwieldy, has
    # been meticulously timed to be the most space and time efficient... Yup!
    for item_first in itermod.iter_items(*iterables):
        break

    # Return this item.
    return item_first

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


@type_check
def get_item_first_satisfying_or_sentinel(
    iterable: IterableTypes, predicate: CallableTypes) -> object:
    '''
    First item of the passed iterable satisfying the passed **predicate**
    (i.e., callable accepting one parameter, returning ``True`` only if this
    parameter suffices) if this iterable contains such an item *or* the
    **sentinel singleton** (i.e.,
    :attr:`betse.util.type.obj.sentinels.SENTINEL`) otherwise.

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
        Either:

        * If one or more items of this iterable satisfy this predicate, the
          first such item.
        * Else, the sentinel singleton.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.obj.sentinels import SENTINEL

    # Collective efficiency is our middle names.
    return next((item for item in iterable if predicate(item)), SENTINEL)

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

# ....................{ GETTERS ~ str                     }....................
#FIXME: Exercise both the "item_key" and "is_defaultable" parameters with
#exhaustive unit tests. This function has become astonishingly complex.
@type_check
def get_item_str_uniquified(
    # Mandatory parameters.
    iterable: IterableTypes,
    item_str_format: str,

    # Mutually exclusive parameters.
    item_attr_name: StrOrNoneTypes = None,
    item_key: StrOrNoneTypes = None,
    is_defaultable: bool = False,
) -> str:
    '''
    Create and return a new machine-readable string guaranteed to be unique
    across each attribute with the passed name of each item of the
    passed iterable if the ``item_attr_name`` parameter is non-``None`` *or*
    across each dictionary key with the passed name of each item of the passed
    iterable if the ``item_key`` parameter is non-``None``.

    Design
    ----------
    This function enforces uniqueness of either an object attribute *or*
    dictionary key assumed to be used by the caller as an SQL-like primary key
    uniquely identifying each item of this iterable. Specifically, this
    function:

    * If the ``item_attr_name`` parameter is non-``None``, requires each item
      of this iterable to declare an attribute with the passed name whose value
      is an arbitrary string.
    * If the ``item_key`` parameter is non-``None``, requires each item of this
      iterable to be a dictionary defining a key with the passed name whose
      value is an arbitrary string.

    Exactly one of the mutually exclusive ``item_attr_name`` and ``item_key``
    parameters *must* be passed. If neither *or* both of these parameters are
    passed, an exception is raised.

    This function then synthesizes a string suitable for use by callers as the
    value of that variable or key for a currently non-existing item of this
    iterable assumed to be created by the caller after calling this function.

    Parameters
    ----------
    iterable : IterableTypes
        Iterable to be inspected.
    item_str_format : str
        Format specifier containing a ``{}`` substring (e.g., ``Item ({}).``),
        iteratively interpolated by this function with an arbitrary integer to
        produce the returned string.
    item_attr_name : StrOrNoneTypes
        Name of the attribute declared by *all* items of this list whose string
        values are to be uniquified. Defaults to ``None``, in which case the
        optional ``item_key`` parameter *must* be non-``None``.
    item_key : StrOrNoneTypes
        Key declared by *all* dictionary items of this list whose string values
        are to be uniquified. Defaults to ``None``, in which case the optional
        ``item_attr_name`` parameter *must* be non-``None``.
    is_defaultable : bool
        ``True`` only if this function implicitly accepts **non-compliant
        items** (i.e., items *not* declaring either an attribute with the name
        ``item_attr_name`` if that parameter is non-``None`` *or* the key
        ``item_key`` if that parameter is non-``None``). Specifically, if this
        boolean is:

        * ``True``, this function silently synthesizes a default string value
          for each non-compliant item of this iterable (e.g.,
          ``item_str_format.format(item_index)``, where ``item_index`` is the
          0-based index of that item in this iterable).
        * ``False``, this function raises an exception on visiting the first
          non-compliant item of this iterable.

        Defaults to ``False``.

    Returns
    ----------
    str
        New machine-readable string guaranteed to both match this format *and*
        be unique across all attributes with this name of all items of
        this iterable.

    Raises
    ----------
    BetseException
        If either:

        * The ``item_attr_name`` parameter is non-``None`` *and*, for one or
          more items of this iterable, either:

          * The ``is_defaultable`` parameter is ``False`` *and* this
            item contains no attribute with this name.
          * The value of this attribute in this item is *not* a string.

        * The ``item_key`` parameter is non-``None`` *and*, for one or more
          items of this iterable, either:

          * This item is *not* a dictionary.
          * The ``is_defaultable`` parameter is ``False`` *and* this
            dictionary does *not* contain this key.
          * The value of this key in this dictionary is *not* a string.
    BetseParamException
        If either:

        * Both of the ``item_attr_name`` *and* ``item_key`` parameters are
          passed.
        * Neither the ``item_attr_name`` *nor* ``item_key`` parameters are
          passed.
    BetseStrException
        If the passed format specifier contains no ``{}`` substring.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.iterable.mapping import mappings
    from betse.util.type.obj import objects
    from betse.util.type.text.string import strs

    # Log this formatting.
    logs.log_debug(
        'Uniquifying iterable item variable "%s" with template "%s"...',
        item_attr_name, item_str_format)

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
    strs.die_unless_substr(text=item_str_format, substr='{}')

    # Callable accepting an arbitrary item and 0-based index of this item in
    # this iterable and returning the string value of the desired attribute or
    # dictionary key from this item. The signature of this callable resembles:
    #     def item_attr_getter(item: object, item_index: int) -> str
    item_str_getter = None

    # If uniquifying an attribute of these items...
    if item_attr_name is not None:
        # If uniquifying a dictionary key of these items, raise an exception.
        # There can only be one true item query paradigm.
        if item_key is not None:
            raise BetseParamException(
                '"item_attr_name" and "item_key" parameters both passed.')

        # If this function implicitly accepts non-compliant items, define a
        # lambda retrieving this string attribute from the passed item in a
        # manner defaulting to a string formatted with this item's index.
        if is_defaultable:
            item_str_getter = lambda item, item_index: (
                objects.get_attr_or_default(
                    obj=item,
                    attr_name=item_attr_name,
                    attr_default=item_str_format.format(item_index),
                    attr_type=str,
                ))
        # Else, this function explicitly rejects non-compliant items. In this
        # case, define a lambda retrieving this string attribute from the
        # passed item *WITHOUT* defaulting if this attribute is undefined.
        else:
            item_str_getter = lambda item, item_index: objects.get_attr(
                obj=item,
                attr_name=item_attr_name,
                attr_type=str,
            )
    # Else if uniquifying a dictionary key of these items...
    elif item_key is not None:
        # If this function implicitly accepts non-compliant items, define a
        # lambda retrieving this string value from the passed dictionary in a
        # manner defaulting to a string formatted with this dictionary's index.
        if is_defaultable:
            item_str_getter = lambda item, item_index: (
                mappings.get_key_value_or_default(
                    mapping=item,
                    key=item_key,
                    value_default=item_str_format.format(item_index),
                    value_type=str,
                ))
        # Else, this function explicitly rejects non-compliant items. In this
        # case, define a lambda retrieving this string value from the
        # passed dictionary *WITHOUT* defaulting if this key is undefined.
        else:
            item_str_getter = lambda item, item_index: mappings.get_key_value(
                mapping=item,
                key=item_key,
                value_type=str,
            )
    # Else, neither an attribute nor dictionary key of these items is
    # being uniquified. In this case, raise an exception.
    else:
        raise BetseParamException(
            '"item_attr_name" and "item_key" parameters not passed.')

    # Arbitrary integer to be formatted into the string to be returned. To
    # reduce the likelihood of collisions with existing items of this iterable,
    # this integer defaults to the 0-based index of the next item to add to
    # this iterable.
    #
    # Iteration below increments this integer until the string created by
    # formatting this integer into the passed format specifier produces a
    # string guaranteed to be unique across all items of this iterable.
    item_id = len(iterable) + 1

    # Uniquified string to be created and returned by this function.
    item_str = None

    # "True" only if the "item_str" string still collides with at least one
    # string value of an attribute with this name of an item of this
    # iterable and hence has yet to be uniquified.
    is_item_str_collides = True

    # Unordered set of all string values of each desired attribute or
    # dictionary key  of all items of this iterable if each item declares such
    # a variable or key *OR* raise an exception otherwise.
    item_strs = set(
        item_str_getter(item, item_index)
        for item_index, item in enumerate(iterable))

    # While the string value to be returned still collides with at least one
    # string value of an attribute with this name of an item of this
    # iterable and hence has yet to be uniquified...
    while is_item_str_collides:
        # Uniquified string to be created and returned by this function.
        item_str = item_str_format.format(item_id)

        # If this value collides with the uniquified string to be created
        # and returned by this function, this string is non-unique. In this
        # case, increment this identifier, format a new attribute with this
        # identifier, and continue searching.
        is_item_str_collides = item_str in item_strs

        # Unconditionally increment the integer to be formatted into the
        # "item_str" string by the next iteration of this loop.
        item_id += 1

    # Return this uniquified string.
    return item_str
