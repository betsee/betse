#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level **non-string iterable** (i.e., non-string object implementing the
abstract :class:`collections.abc.Iterable` base class) facilities.

See Also
----------
:func:`betse.util.type.types.is_iterable`
    Further details on what constitutes iterables and non-string iterables.
'''

# ....................{ IMPORTS                           }....................
import itertools, pprint
from betse.exceptions import BetseIterableException
from betse.util.type import types
from betse.util.type.types import (
    type_check,
    ClassOrNoneTypes,
    GeneratorType,
    IterableTypes,
    MappingType,
    SizedType,
)
from collections import deque

# ....................{ CONSUMERS                         }....................
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
        strictly positive (i.e., ``iterations >= 1``).

    Returns
    ----------
    object
        Object yielded by the last iterable iteration (i.e., the last
        ``next()`` method called on this iterable) if any *or* ``None`` if this
        iterable was already exhausted (i.e., empty) when passed.

    See Also
    ----------
    http://docs.python.org/3/library/itertools.html
        Official documentation strongly inspiring this function.
    '''

    #FIXME: Convert into a propert exception.
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
    * Does *not* explicitly halt with a :class:`StopIteration` exception and
      hence is infinite, this function reduces to an **infinite loop.**

    Parameters
    ----------
    iterable : IterableTypes
        Iterable to be exhausted.

    Returns
    ----------
    object
        Object yielded by the last iterable iteration (i.e., the last
        ``next()`` method called on this iterable) if any *or* ``None`` if this
        iterable was already exhausted (i.e., empty) when passed.

    See Also
    ----------
    http://docs.python.org/3/library/itertools.html
        Official documentation strongly inspiring this function.
    '''

    # For efficiency, feed this iterable into a zero-length deque retaining
    # only the last iterated value if any.
    iterable_deque = deque(iterable, maxlen=1)

    # If this iterable was *NOT* already exhausted, return its last value.
    if iterable_deque:
        return iterable_deque[0]
    # Else, this iterable was already exhausted. Return nothing.
    else:
        return None

# ....................{ CONVERTERS                        }....................
@type_check
def to_iterable(
    # Mandatory parameters.
    iterable: IterableTypes,

    # Optional parameters.
    cls: ClassOrNoneTypes = None,
    item_cls: ClassOrNoneTypes = None,
) -> IterableTypes:
    '''
    Convert the passed input iterable into an output iterable of the passed
    iterable type and/or convert each item of this input iterable into an item
    of this output iterable of the passed item type.

    Specifically, if this iterable type is:

    * Either ``None`` *or* of the same type as the passed iterable type, then:

      * If this item type is non-``None``, this input iterable is converted
        into an instance of the same iterable type whose items are converted
        into instances of this item type.
      * Else, this input iterable is returned as is (i.e., unmodified).

    * A non-Numpy iterable (e.g., :class:`list`) and the passed iterable type
      is that of:

      * Another non-Numpy iterable (e.g., :class:`tuple`), then:

        * If this item type is non-``None``, this input iterable is converted
          into an instance of this iterable type whose items are converted into
          instances of this item type.
        * Else, this input iterable is merely converted into an instance of
          this iterable type.

        In either case, this output iterable's ``__init__`` method is required
        to accept this input iterable as a single positional argument.

      * A Numpy array, then:

        * If this item type is non-``None``, an exception is raised. Numpy
          scalar types map poorly to Python scalar types.
        * Else, this input iterable is converted to a Numpy array via the
          :func:`betse.lib.numpy.nparray.from_iterable` function.

    * A Numpy array, then:

      * If this item type is non-``None``, an exception is raised. Numpy
        scalar types map poorly to Python scalar types.
      * Else, this input iterable is converted to a Numpy array via the
        :func:`betse.lib.numpy.nparray.from_iterable` function.

    Parameters
    ----------
    iterable: IterableTypes
        Input iterable to be converted.
    cls : ClassOrNoneTypes
        Type of the output iterable to convert this input iterable into.
        Defaults to ``None``, in which case the same type as that of this input
        iterable is defaulted to.
    item_cls : ClassOrNoneTypes
        Type to convert each item of the output iterable into. Defaults to
        ``None``, in which case these items are preserved as is.

    Returns
    ----------
    IterableTypes
        Output iterable converted from this input iterable.

    Raises
    ----------
    BetseIterableException
        If converting either to or from a Numpy array *and* this item type is
        non-``None``.
    '''

    # Avoid importing third-party packages at the top level, for safety.
    from betse.lib.numpy import nparray
    from betse.util.type.cls import classes
    from betse.util.type.iterable import generators
    from numpy import ndarray

    # Type of the input iterable.
    iterable_src_cls = type(iterable)

    # Type of the output iterable.
    iterable_trg_cls = cls

    # If the caller requested no explicit type conversion...
    if iterable_trg_cls is None:
        # If the input iterable is a generator, default the type of the output
        # iterable to the optimally space- and time-efficient iterable: tuple.
        # Why? Because generators *CANNOT* be explicitly instantiated.
        if generators.is_generator(iterable):
            iterable_trg_cls = tuple
        # Else, the input iterable is *NOT* a generator and hence is assumed to
        # be of a standard type that *CAN* be explicitly instantiated. In this
        # case, default the type of the output iterable to this type.
        else:
            iterable_trg_cls = iterable_src_cls

    # if cls is not None else iterable_src_cls

    # If the input and output iterables are of the same type *AND* no item
    # conversion was requested, efficiently reduce to a noop.
    if iterable_src_cls is iterable_trg_cls and item_cls is None:
        return iterable

    # Else if the input iterable is a Numpy array...
    if iterable_src_cls is ndarray:
        # If the caller requested item conversion, raise an exception.
        if item_cls is not None:
            raise BetseIterableException(
                'Numpy array not convertible to item class "{}".'.format(
                    classes.get_name_unqualified(item_cls)))

        # Defer to logic elsewhere.
        return nparray.to_iterable(array=iterable, cls=iterable_trg_cls)

    # Else if the output iterable is a Numpy array, defer to logic elsewhere.
    if iterable_trg_cls is ndarray:
        # If the caller requested item conversion, raise an exception.
        if item_cls is not None:
            raise BetseIterableException(
                'Numpy array not convertible to item class "{}".'.format(
                    classes.get_name_unqualified(item_cls)))

        return nparray.from_iterable(iterable)

    # Else, neither the input or output iterables are Numpy arrays. In this
    # case, return an output iterable of the desired type containing the items
    # of this input iterable either...
    return (
        # Unmodified if no item conversion was requested *OR*...
        iterable_trg_cls(iterable) if item_cls is None else
        # Converted to this item type otherwise.
        iterable_trg_cls(item_cls(item) for item in iterable))


@type_check
def to_str(iterable: IterableTypes) -> str:
    '''
    Convert the passed input iterable into an output string intended to be
    substantially more human-readable than the standard implementations of
    most ``__str__`` special methods implicitly invoked by the
    :meth:`str.__init__` method.

    Parameters
    ----------
    iterable: IterableTypes
        Input iterable to be converted.

    Returns
    ----------
    str
        Output string converted from this input iterable.
    '''

    # Simplicity is a place in this codebase.
    return pprint.pformat(iterable)

# ....................{ INVERTERS                         }....................
@type_check
def invert_iterable_unique(iterable: IterableTypes) -> MappingType:
    '''
    Dictionary inverted from the passed iterable if **internally unique**
    (i.e., containing no duplicate items) *or* raise an exception otherwise.

    Specifically:

    * If this iterable is a dictionary, the
      :func:`betse.util.type.iterable.mapping.mappings.invert_dict_unique` function is
      silently deferred to. The type of the returned dictionary is guaranteed
      to be the same as the type of the passed dictionary.
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
    from betse.util.type.iterable import itertest
    from betse.util.type.iterable.mapping import mappings

    # If this iterable is a mapping...
    if mappings.is_mapping(iterable):
        # Invert this iterable specifically as a mapping.
        #
        # While mappings are technically iterables and hence invertable via the
        # generic approach applied below, doing so incorrectly returns a
        # dictionary mapping from the keys of this dictionary to arbitrary
        # 0-based integers. As the adjective "arbitrary" implies, this renders
        # the resulting dictionary effectively useless for most purposes.
        return mappings.invert_map_unique(iterable)
    # Else, this iterable is a non-mapping. In this case, a generic approach
    # suffices... usually.
    else:
        # If one or more items of this iterable are duplicates, raise an
        # exception.
        itertest.die_unless_items_unique(iterable)

        # One-liners for Great Glory.
        return {item: item_index for item_index, item in enumerate(iterable)}

# ....................{ ITERATORS                         }....................
@type_check
def iter_items(*iterables: IterableTypes) -> GeneratorType:
    '''
    Generator yielding each item of each of the passed iterables (in both the
    internal order of each iterable and the passed order of iterables),
    effectively "chaining" these iterables together.

    This function avoids instantiating any objects (other than the returned
    generator) and hence is *probably* the maximally space and time efficient
    implementation of this task. In particular, this function is preferable for
    pure iteration over multiple iterables, in which case a composite iterable
    of the same or different type is neither required nor desired.

    Parameters
    ----------
    iterables : tuple[IterableTypes]
        Tuple of all iterables to be iterated over.

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
        >>> from betse.util.type.iterable import iterables
        >>> in_xanadu = ('did', 'Kubla', 'Khan')
        >>> a_stately = ('pleasure-dome', 'decree:')
        >>> for word in iterables.iter_items(in_xanadu, a_stately):
        ...     print(word, end=' ')
        did Kubla Khan pleasure-dome decree:
    '''

    # Pure Python 3.x. Join the party.
    #
    # Note that the itertools.chain.from_iterable() could also be called here,
    # but that doing so would be both less Pythonic *AND* no more efficient
    # than the current approach. Ergo, simplicity wins.
    for iterable in iterables:
        yield from iterable

# ....................{ JOINERS                           }....................
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
        >>> from betse.util.type.iterable import iterables
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

# ....................{ REVERSERS                         }....................
@type_check
def reverse(iterable: IterableTypes) -> IterableTypes:
    '''
    New iterable of possibly differing type guaranteed to be the **reverse**
    (i.e., contain all items in reverse order) of the passed iterable.

    The type of the iterable created and returned by this function is either:

    * If the passed iterable is *not* reversible as is (i.e., the
      :func:`is_reversible` function returns ``False`` for this iterable), a
      :class:`tuple` -- the most space- and time-efficient container capable of
      containing *any* items.
    * If the passed iterable is reversible as is (i.e., the
      :func:`is_reversible` function returns ``True`` for this iterable), the
      same type as this iterable (i.e., ``type(iterable)``).

    Caveats
    ----------
    **This high-level reversal function should typically be called in lieu of
    the low-lovel :func:`reversed` builtin.** While this function generically
    supports *all* possible iterable types, :func:`reversed` only specifically
    supports **reversible iterable types** (i.e., instances of which the
    :func:`is_reversible` function returns ``True`` for).

    **The new iterable created and returned by this function is only a shallow
    copy of the passed iterable for efficiency.** These two iterables thus
    share the same items in the reverse order. Callers instead requiring a
    deep copy must do so either before or after calling this function.

    Parameters
    ----------
    iterable : IterableTypes
        Iterable to be reversed. For generality, this iterable remains
        unchanged by this function.

    Returns
    ----------
    IterableTypes
        Iterable reversed from the passed iterable.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.iterable import itertest

    # If this iterable is *NOT* reversible as is, convert this iterable into
    # the most space- and time-efficient iterable containing the same items
    # that *IS* reversible -- in this case, a tuple.
    if not itertest.is_reversible(iterable):
        iterable = tuple(iterable)

    # Return the result of the efficient reversed() builtin on this iterable,
    # now guaranteed to be reversible as is.
    return reversed(iterable)

# ....................{ ZIPPERS                           }....................
#FIXME: Unit test us up.
@type_check
def zip_isometric(*iterables: IterableTypes) -> GeneratorType:
    '''
    Generator zipping all passed iterables required to be of the same length.

    This generator iteratively yields an n-tuple, where:

    * n is the length of each passed iterable.
    * The i-th element of this tuple is in the i-th passed iterable.

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
        If any passed iterable differs in length from any other passed
        iterable.

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
    iterables identified by the passed `n`-tuple is of smaller length than
    other iterables in this tuple of iterables.

    This private function is *only* intended to be called by the
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

    # If this erroneously short iterable is of predefined length (e.g., is
    # *NOT* a generator of dynamic length), begin this exception message with
    # this length and end this message with this iterable's contents.
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
