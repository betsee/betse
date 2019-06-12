#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level **mapping utilities** (i.e., functions operating on dictionary-like
types and instances).
'''

# ....................{ IMPORTS                           }....................
from betse.exceptions import BetseMappingException, BetseMappingKeyException
from betse.util.type.types import (
    type_check,
    MappingType,
    HashableType,
    IterableTypes,
    TestableOrNoneTypes,
)
from copy import deepcopy

# ....................{ EXCEPTIONS                        }....................
@type_check
def die_unless_keys_equal(*mappings: MappingType) -> None:
    '''
    Raise an exception unless all of the passed dictionaries contain the exact
    same keys.

    Equivalently, this function raises an exception if any key of any passed
    dictionary is *not* a key of any other such dictionary.

    Parameters
    ----------
    mapping : MappingType
        Dictionary to be inspected.

    Raises
    ----------
    BetseMappingException
        If any key of any passed dictionary is *not* a key of any other such
        dictionary.

    See Also
    ----------
    :func:`is_keys_equal`
        Further details.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.text.string import strjoin

    # If one or more of these dictionaries contain differing keys...
    if not is_keys_equal(*mappings):
        # First passed mapping. Since the is_keys_equal() function necessarily
        # returns true if either no mappings or only one mapping are passed,
        # this function returning false implies that two or more mappings are
        # passed. Ergo, this mapping is guaranteed to exist.
        mapping_first = mappings[0]

        # For each mapping excluding the first...
        for mapping in mappings[1:]:
            # If the keys of this mapping differ from those of the first...
            if not is_keys_equal(mapping, mapping_first):
                # Set of all keys differing between these two mappings.
                keys_unequal = mapping.keys().symmetric_difference(
                    mapping_first.keys())

                # Grammatically correct noun describing the number of such
                # keys.
                keys_noun = 'key' if len(keys_unequal) == 1 else 'keys'

                # Raise an exception embedding this set.
                raise BetseMappingException(
                    'Dictionary {} {} differ.'.format(
                        keys_noun,
                        strjoin.join_as_conjunction_double_quoted(
                            *keys_unequal)))


@type_check
def die_unless_values_unique(mapping: MappingType) -> None:
    '''
    Raise an exception unless all values of the passed dictionary are unique.

    Equivalently, this function raises an exception if any two key-value pairs
    of this dictionary share the same values.

    Parameters
    ----------
    mapping : MappingType
        Dictionary to be inspected.

    Raises
    ----------
    BetseMappingException
        If at least one value of this dictionary is a duplicate.

    See Also
    ----------
    :func:`is_values_unique`
        Further details.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.iterable import iterget
    from betse.util.type.text.string import strjoin

    # If one or more values of this dictionary are duplicates...
    if not is_values_unique(mapping):
        # Set of all duplicate values in this dictionary.
        values_duplicate = iterget.get_items_duplicate(mapping.values())

        # Grammatically correct noun describing the number of such values.
        values_noun = 'value' if len(values_duplicate) == 1 else 'values'

        # Raise an exception embedding this set.
        raise BetseMappingException(
            'Dictionary {} {} duplicate.'.format(
                values_noun,
                strjoin.join_as_conjunction_double_quoted(*values_duplicate)))

# ....................{ TESTERS ~ key                     }....................
@type_check
def is_key(mapping: MappingType, *keys: HashableType) -> bool:
    '''
    ``True`` only if the passed dictionary contains *all* passed keys.

    Parameters
    ----------
    mapping : MappingType
        Dictionary to be tested.
    keys : tuple[HashableType]
        Tuple of all keys to be tested for.

    Returns
    ----------
    bool
        ``True`` only if this dictionary contains *all* passed keys.
    '''

    return (
        # If only one key is passed, optimize this common edge case with the
        # standard idiom for testing key containment.
        keys[0] in mapping if len(keys) == 1 else
        # Else, two or more keys are passed. In this case, fallback to a
        # general-case strategy testing key containment in a single one-liner.
        # And yes: this is ridiculously awesome.
        set(keys).issubset(mapping)
    )


@type_check
def is_keys_equal(*mappings: MappingType) -> bool:
    '''
    ``True`` only if all passed dictionaries contain the exact same keys.

    Parameters
    ----------
    mappings : Tuple[MappingType]
        Tuple of all dictionaries to be tested.

    Returns
    ----------
    bool
        ``True`` only if these dictionaries contain the exact same keys.
    '''

    # If two mappings are passed, prematurely optimize this common case by
    # directly testing the keys contained in these two mappings for equality.
    if len(mappings) == 2:
        return mappings[0].keys() == mappings[1].keys()
    # Else if either no mappings or only one mapping are passed, return true.
    # Why? Because:
    #
    # * If no mappings are passed, this edge case is functionally equivalent to
    #   the edge case in which one or more empty mappings (i.e., mappings
    #   containing no key-value pairs) are passed. Since the key containers for
    #   empty mappings are themselves empty, these key containers contain the
    #   exact same keys -- namely, none.
    # * If only one mapping is passed, this mapping by definition contains the
    #   exact same keys as itself.
    elif len(mappings) < 2:
        return True
    # Else, three or more mappings are passed. In this case, defer to a general
    # case algorithm.

    # Keys of the first passed mapping.
    mapping_first_keys = mappings[0].keys()

    # Return true only if...
    return all(
        # This mapping contains the same keys as the first such mapping.
        mapping.keys() == mapping_first_keys
        # For each mapping excluding the first.
        for mapping in mappings[1:]
    )

# ....................{ TESTERS ~ value                   }....................
@type_check
def is_values_unique(mapping: MappingType) -> bool:
    '''
    ``True`` only if all values of the passed dictionary are **unique** (i.e.,
    if *no* two key-value pairs of this dictionary share the same values).

    Parameters
    ----------
    mapping : MappingType
        Dictionary to be inspected.

    Returns
    ----------
    bool
        ``True`` only if *all* values of this dictionary are unique.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.iterable import itertest

    # For sanity, defer to an existing low-level tester.
    return itertest.is_items_unique(mapping.values())

# ....................{ GETTERS                           }....................
@type_check
def get_key_value(
    mapping: MappingType, key: HashableType, **kwargs) -> object:
    '''
    Value of the passed key in the passed mapping if this mapping contains this
    key *or* raise an exception otherwise (i.e., if this mapping contains no
    such key), optionally validated to be of the passed type.

    Parameters
    ----------
    mapping : MappingType
        Dictionary to be inspected.
    key : HashableType
        Key whose value is to be retrieved.

    All remaining keyword arguments are passed as is to the
    :func:`get_key_value_or_default` function.

    Returns
    ----------
    object
        Value of this key in this mapping.

    Raises
    ----------
    BetseMappingKeyException
        If this mapping contains no such key.
    BetseTypeException
        If the ``value_type`` parameter is non-``None`` and the type of the
        current value of this key is *not* an instance of ``value_type``.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.obj.sentinels import SENTINEL

    # Value of this key in this mapping if any *OR* the sentinel otherwise.
    key_value = get_key_value_or_sentinel(
        mapping=mapping, key=key, **kwargs)

    # If this mapping contains no such key, raise an exception.
    if key_value is SENTINEL:
        raise BetseMappingKeyException(
            'Mapping key "{}" not found.'.format(key))
    # Else, this mapping contains this key.

    # Return this value.
    return key_value


@type_check
def get_key_value_or_sentinel(
    mapping: MappingType, key: HashableType, **kwargs) -> object:
    '''
    Value of the passed key in the passed mapping if this mapping contains this
    key *or* the sentinel singleton otherwise (i.e., if this mapping contains
    no such key), optionally validated to be of the passed type.

    This function enables callers to safely distinguish between non-existing
    keys and existing keys whose values are ``None``.

    Parameters
    ----------
    mapping : MappingType
        Dictionary to be inspected.
    key : HashableType
        Key whose value is to be retrieved.

    All remaining keyword arguments are passed as is to the
    :func:`get_key_value_or_default` function.

    Returns
    ----------
    object
        Either:

        * If this dictionary contains this key, this key's value.
        * Else, the **sentinel singleton** (i.e.,
          :attr:`betse.util.type.obj.sentinels.SENTINEL`).

    Raises
    ----------
    BetseTypeException
        If the ``value_type`` parameter is non-``None`` and the type of the
        current value of this key is *not* an instance of ``value_type``.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.obj.sentinels import SENTINEL

    # Return the current value of this key in this mapping if any *OR* the
    # sentinel otherwise.
    return get_key_value_or_default(
        mapping=mapping, key=key, value_default=SENTINEL, **kwargs)


@type_check
def get_key_value_or_default(
    # Mandatory parameters.
    mapping: MappingType,
    key: HashableType,
    value_default: object,

    # Optional parameters.
    value_type: TestableOrNoneTypes = None,
) -> object:
    '''
    Value of the passed key in the passed mapping if this mapping contains this
    key *or* the passed default value otherwise (i.e., if this mapping contains
    no such key), optionally validated to be of the passed type.

    Parameters
    ----------
    mapping : MappingType
        Dictionary to be inspected.
    key : HashableType
        Key to return the current value of.
    value_default : object
        Default value to be returned if this dictionary contains no such key.
    value_type : TestableOrNoneTypes
        Expected type of the current value of this key. This function
        effectively performs the equivalent of the :meth:`type_check` decorator
        at runtime by raising an exception if all of the following apply:

        * This type is *not* ``None``.
        * This value is *not* this default value, implying this dictionary to
          contain this key.
        * This value is *not* an instance of this type.

        Defaults to ``None``, in which case no such type checking is performed.

    Returns
    ----------
    object
        Either:

        * If this dictionary contains this key, this key's value.
        * Else, this default value.

    Raises
    ----------
    BetseTypeException
        If the ``value_type`` parameter is non-``None`` and the type of the
        current value of this key is *not* an instance of ``value_type``.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.obj import objtest

    # Value of this key in this mapping if any *OR* this default value.
    key_value = mapping.get(key, value_default)

    # If this value is to be type-checked *AND* is *NOT* this default value
    # (which by definition already satisfies caller requirements regardless of
    # type), type-check this value.
    if value_type is not None and key_value is not value_default:
        objtest.die_unless_instance(obj=key_value, cls=value_type)

    # Return this value.
    return key_value

# ....................{ COPIERS                           }....................
@type_check
def copy_deep(mapping: MappingType) -> MappingType:
    '''
    Dictionary of all key-value pairs deeply (i.e., recursively) duplicated
    from the passed dictionary.

    This function should *always* be called in lieu of the standard
    :meth:`dict.__init__` and :meth:`dict.copy` methods, which only perform
    shallow dictionary copies. These copies fail to copy data structures nested
    in the values of the original dictionary, inviting subtle synchronization
    woes on subsequently modifying either the original or copied dictionaries.

    Parameters
    ----------
    mapping: MappingType
        Dictionary to be deeply copied.

    Returns
    ----------
    MappingType
        Dictionary of all key-value pairs deeply (i.e., recursively) duplicated
        from the passed dictionary.
    '''

    #FIXME: Does this simplistic approach guarantee the returned mapping to be
    #of the same type as the passed mapping?
    return deepcopy(mapping)


#FIXME: Well, this is rather awkward. Rather than define a completely separate
#function, it would be dramatically preferable to simply pass a new optional
#"keys_remove: IterableOrNoneTypes = None" parameter to the copy_deep()
#function defined above. When this parameter is:
#
#* "None", the copy_deep() function should reduce to its current one-liner.
#* Non-"None", the copy_deep() function should generalize to this function's
#  current implementation -- with the obvious caveat that the remove_key()
#  function should be generalized to accept an iterable of keys to be removed.
@type_check
def copy_map_sans_key(mapping: MappingType, key: HashableType) -> MappingType:
    '''
    Dictionary of all key-value pairs excluding that whose key is the passed
    key deeply (i.e., recursively) duplicated from the passed dictionary.

    Parameters
    ----------
    mapping: MappingType
        Dictionary to be deeply copied.
    key : HashableType
        Key to be removed from this dictionary.

    Returns
    ----------
    MappingType
        Dictionary of all key-value pairs excluding that whose key is this
        key deeply (i.e., recursively) duplicated from this dictionary.

    Raises
    ----------
    :class:`KeyError`
        If this dictionary contains no such key.

    See Also
    ----------
    :func:`copy_deep`
        Further details on map copying.
    :func:`remove_key`
        Further details on key removal.
    '''

    # Deep copy of this dictionary.
    mapping_copy = copy_deep(mapping)

    # Remove this key from this copy in-place.
    remove_key(mapping=mapping_copy, key=key)

    # Return this copy.
    return mapping_copy

# ....................{ INVERTERS                         }....................
#FIXME: Rename to simply invert_unique().
@type_check
def invert_map_unique(mapping: MappingType) -> MappingType:
    '''
    Dictionary inverted from the passed dictionary if no two key-value pairs of
    this dictionary share the same values *or* raise an exception otherwise.

    Specifically, the returned dictionary maps from each value to each key of
    the passed dictionary *and* is guaranteed to be the same type as that of
    the passed dictionary.

    Parameters
    ----------
    mapping : MappingType
        Dictionary to be inverted. The type of this dictionary *must* define an
        ``__init__`` method accepting a single parameter whose value is an
        iterable of 2-iterables ``(key, value)`` providing all key-value pairs
        with which to initialize a new such dictionary. See the
        :meth:`dict.__init__` method for further details.

    Returns
    ----------
    MappingType
        Dictionary inverted from this dictionary as detailed above.

    Raises
    ----------
    BetseMappingException
        If one or more key-value pairs of this dictionary share the same
        values.

    See Also
    ----------
    https://stackoverflow.com/a/1679702/2809027
        StackOverflow answer strongly inspiring this implementation.
    '''

    # If any values of this dictionary are are duplicates, raise an exception.
    die_unless_values_unique(mapping)

    # Type of this dictionary.
    mapping_type = type(mapping)

    # If this is an unordered dictionary, return a dictionary comprehension
    # efficiently inverting this dictionary in the trivial way.
    if mapping_type is dict:
        return {value: key for key, value in mapping.items()}
    # Else, this is possibly an ordered dictionary. In this case, a
    # considerably less trivial and slightly less efficient approach is needed.
    else:
        # Iterable of reversed 2-iterables "(value, pair)" for each key-value
        # pair of the passed dictionary. Dismantled, this is:
        #
        # * "mapping.items()", an iterable of 2-iterables "(key, value)" for
        #   each key-value pair of the passed dictionary.
        # * "reversed", a builtin which when passed such a 2-iterable returns
        #   the reversed 2-iterable "(value, pair)" for that key-value pair.
        # * "map(...)", a builtin applying the prior builtin to each such pair.
        value_key_pairs = map(reversed, mapping.items())

        # Return a new instance of this type of dictionary by invoking the
        # "dict(iterable)" form of this type's __init__() method. To quote the
        # dict.__init__() docstring:
        #
        # "dict(iterable) -> new dictionary initialized as if via:
        #      d = {}
        #      for k, v in iterable:
        #          d[k] = v"
        return mapping_type(value_key_pairs)

# ....................{ MERGERS                           }....................
#FIXME: Refactor as follows:
#
#* Shift into a new "mapmerge" submodule of the same subpackage.
#* Define a new "MergeCollisionType" enumeration in this submodule supporting
#  at least the following three members:
#  * "RAISE_EXCEPTION".
#  * "PREFER_FIRST", giving higher precedence to keys in mappings passed
#    earlier.
#  * "PREFER_LAST", giving higher precedence to keys in mappings passed later.
#* Add a new optional parameter with the following signature:
#    on_collision: MergeCollisionType = MergeCollisionType.RAISE_EXCEPTION,
#* Implement these collision strategies in this function. We've already
#  implemented the "PREFER_LAST" strategy below. The "PREFER_FIRST" strategy
#  probably just reduces to iteritavely calling the dict.update() method on the
#  newly created dictionary to be returned in a reasonably intelligent manner.
#  The "RAISE_EXCEPTION" strategy may require the most work, but should prove
#  mostly trivial.
@type_check
def merge_maps(
    # Mandatory parameters.
    mappings: IterableTypes,

    # Optional parameters.
) -> MappingType:
    '''
    Dictionary of all key-value pairs deeply (i.e., recursively) merged
    together from all passed dictionaries (in the passed order).

    Caveats
    ----------
    If the ``on_collision`` parameter is:

    * :attr:`MergeCollisionType.`, **order is insignificant.** In this case, this
      function raises an exception if any two of the passed dictionaries
      **collide** (i.e., define the same key). Since this prevents key
      collisions, *no* implicit precedence exists between these dictionaries.
    * Either :attr:`MergeCollisionType.`, **order is significant.** In this case, dictionaries passed later take precedence over
      dictionaries passed earlier. Ergo, the last passed dictionary takes
      precedence over *all* other passed dictionaries. Whenever any two passed
      dictionaries collide (i.e., contain the same key), the returned dictionary
      contains a key-value pair for that key whose value is that of the key-value
      pair for the same key of whichever of the two dictionaries was passed last.

    Parameters
    ----------
    mappings : Tuple[MappingType]
        Tuple of all dictionaries to be merged.

    Returns
    ----------
    MappingType
        Dictionary merged from and of the same type as the passed dictionaries.
        Note lastly that the class of the passed dictionary *must* define an
        ``__init__()`` method accepting a dictionary comprehension.

    See Also
    ----------
    :meth:`dict.update`
        Standard method merging two dictionaries, which should typically be
        called instead of this slower function in this specific use case.
    http://treyhunner.com/2016/02/how-to-merge-dictionaries-in-python
        Blog post strongly inspiring this implementation. Thanks, Trey!
    '''

    # Avoid circular import dependencies.
    from betse.util.type.iterable import itertest, sequences

    # If no mappings were passed, raise an exception.
    sequences.die_if_empty(mappings, label='Mapping')

    # If any of the passed mappings is *NOT* a mapping, raise an exception.
    itertest.die_unless_items_instance_of(iterable=mappings, cls=MappingType)

    # Type of dictionary to be returned.
    dict_type = type(mappings[0])

    # Dictionary merged from the passed dictionaries via a doubly-nested
    # dictionary comprehension. While there exist a countably infinite number
    # of approaches to merging dictionaries in Python, this approach is known
    # to be the most efficient for general-purpose merging of arbitrarily many
    # dictionaries under Python >= 3.4. See also Trey Hunter's exhaustive
    # commentary (complete with timings) at the above URL.
    dict_merged = {
        # For safety, deeply copy rather than reuse this value.
        key: deepcopy(value)
        for mapping in mappings
        for key, value in mapping.items()
    }

    # Return a dictionary of this type converted from this dictionary. If the
    # desired type is a "dict", this dictionary is returned as is; else, this
    # dictionary is converted into an instance of the desired type.
    return dict_merged if dict_type is dict else dict_type(dict_merged)

# ....................{ REMOVERS                          }....................
@type_check
def remove_key(mapping: MappingType, key: HashableType) -> None:
    '''
    Remove the key-value pair whose key is the passed key from the passed
    dictionary **in-place** (i.e., by modifying this dictionary rather than
    creating and returning a new dictionary with this key removed) if this
    dictionary contains this key *or* raise an exception otherwise.

    This function is a caller convenience improving codebase readability and
    efficiency. Although there exist multiple means of removing key-value pairs
    from dictionaries, this function implements the most efficient approach.
    These include:

    * The ``del mapping[key]`` idiom, known to be the most efficient approach.
    * The :meth:`dict.pop` method, known to be slightly less efficient than the
      idiomatic approach.

    Parameters
    ----------
    mapping : MappingType
        Dictionary to remove this key from.
    key : HashableType
        Key to be removed from this dictionary.

    Raises
    ----------
    :class:`KeyError`
        If this dictionary contains no such key.

    See Also
    ----------
    :func:`copy_map_sans_key`
        Function creating and returning a new dictionary with this key removed.
    '''

    # The best things in life are free.
    del mapping[key]

# ....................{ ITERATORS                         }....................
@type_check
def iter_keys_ascending(mapping: MappingType) -> list:
    '''
    Iterable of all keys of the passed dictionary sorted in ascending order.

    Parameters
    ----------
    mapping : MappingType
        Dictionary whose keys are to be sorted.

    Returns
    ----------
    list
        Iterable of all keys of this dictionary sorted in ascending order.

    See Also
    ----------
    :func:`betse.util.type.iterable.itersort.sort_ascending`
        Further details on ascending sorting.
    '''

    # Return a list of all keys of this dictionary sorted in ascending order.
    #
    # Note that the itersort.sort_ascending() function is intentionally *NOT*
    # called here, as that function currently fails to support dictionary views
    # in an efficient manner. Naturally, this should be rectified at some time.
    return sorted(mapping)
