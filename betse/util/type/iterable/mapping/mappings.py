#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level **mapping** (i.e., :class:`dict`-like types or instances)
functionality.
'''

# ....................{ IMPORTS                           }....................
from betse.exceptions import BetseMappingKeyException
from betse.util.type.types import (
    type_check, MappingType, HashableType, TestableOrNoneTypes)
from copy import deepcopy

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

# ....................{ CONVERTERS                        }....................
@type_check
def to_str_flat(mapping: MappingType) -> str:
    '''
    Human-readable string flattening all key-value pairs of the passed
    dictionary.

    For readability, each such pair is sorted in ascending lexicographic order
    and formatted as ``{key}: {value}\n``. Ergo:

    * Each key and value is reduced to a string.
    * Each key is delimited from its value by a colon.
    * Each key-value pair is delimited by a newline.

    Caveats
    ----------
    **This function assumes all keys and values of this dictionary to be
    trivially reducible to terse strings,** where "terse" typically implies
    these strings to be no longer than a standard terminal width of 80
    characters. When this is *not* the case, the returned string is likely to
    be non-human-readable.

    Parameters
    ----------
    mapping: MappingType
        Dictionary to be flattened.

    Returns
    ----------
    str
        Human-readable string flattening all key-value pairs of this
        dictionary.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.text.string import strjoin

    # Return the human-readable string produced by joining on newline a
    # generator comprehension yielding the colon-delimited name and value of
    # all environment variables sorted in lexicographic order.
    #
    # Note that the "environ" object is of non-standard type "os._Environ",
    # which the pprint.pformat() fails to recognize and hence format as a
    # "dict"-compatible mapping. Ergo, passing "environ" directly to the
    # iterables.to_str() function would produce a non-human-readable string.
    # While this can, of course, be ameliorated by converting "environ" to a
    # "dict" first (e.g., "iterables.to_str(dict(environ))"), doing so still
    # produces less human-readable output than the current approach.
    # return iterables.to_str(dict(environ))
    return strjoin.join_on_newline(
        '{}: {}'.format(key, mapping[key])
        for key in iter_keys_ascending(mapping)
    )

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

    # Avoid circular import dependencies.
    from betse.util.type.iterable.mapping import maptest

    # If any values of this dictionary are are duplicates, raise an exception.
    maptest.die_unless_values_unique(mapping)

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
