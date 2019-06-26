#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level **mapping merger** (i.e., :class:`set`-style union and intersection
operations on dictionary-like types or instances) functionality.
'''

# ....................{ IMPORTS                           }....................
from betse.util.type.enums import make_enum
from betse.util.type.types import type_check, MappingType, IterableTypes
from copy import deepcopy

# ....................{ ENUMERATIONS                      }....................
MergeCollisionPolicy = make_enum(
    class_name='MergeCollisionPolicy',
    member_names=(
        'RAISE_EXCEPTION',
        'PREFER_FIRST',
        'PREFER_LAST',
    ),
)
'''
Enumeration of all supported types of **key collision merger policies** (i.e.,
strategies for merging keys shared between two or more mappings).

Key collisions occur when two or more mappings to be merged contain key-value
pairs containing the same keys but differing values -- or, equivalently, when a
key of a key-value pair of one mapping is also a key of a key-value pair of
another mapping such that the two values of these two pairs differ.

Attributes
----------
RAISE_EXCEPTION : enum
    This policy raises a fatal exception on the first key collision. This
    constitutes the strictest and hence safest such policy.
PREFER_FIRST : enum
    This policy accepts *only* the first key-value pair with this key in these
    mappings and thus ignores all subsequent such pairs in subsequent mappings
    on every key collision. This policy gives higher precedence to keys in
    mappings passed earlier -- the converse of the :attr:`PREFER_LAST` policy.
PREFER_LAST : enum
    This policy accepts *only* the last key-value pair with this key in these
    mappings and thus ignores all prior such pairs in prior mappings on every
    key collision. This policy gives higher precedence to keys in mappings
    passed later -- the converse of the :attr:`PREFER_FIRST` policy.
'''

# ....................{ MERGERS                           }....................
@type_check
def merge_maps(
    # Mandatory parameters.
    mappings: IterableTypes,

    # Optional parameters.
    collision_policy: MergeCollisionPolicy = (
        MergeCollisionPolicy.RAISE_EXCEPTION),
    is_values_copied: bool = False,
) -> MappingType:
    '''
    Dictionary of all key-value pairs merged together from all passed
    dictionaries (in the passed order) with the passed key collision policy.

    Caveats
    ----------
    All classes of all passed dictionaries *must* define ``__init__()`` methods
    accepting a dictionary comprehension.

    Whether the order in which dictionaries are passed is significant
    conditionally depends on the passed ``collision_policy`` parameter.
    Specifically, if this parameter is:

    * :attr:`MergeCollisionPolicy.RAISE_EXCEPTION`, **order is insignificant.**
      In this case, this function simply raises an exception if any two of the
      passed dictionaries **key-collide** (i.e., define key-value pairs of the
      same key but differing values).
    * :attr:`MergeCollisionPolicy.PREFER_FIRST`, **order is significant.** When
      one or more mappings to be merged contain the same key, this function
      accepts only the first key-value pair with this key in these mappings and
      thus ignores all subsequent such pairs in subsequent mappings.
    * :attr:`MergeCollisionPolicy.PREFER_LAST`, **order is significant.** When
      one or more mappings to be merged contain the same key, this function
      accepts only the last key-value pair with this key in these mappings and
      thus ignores all prior such pairs in prior mappings.

    Parameters
    ----------
    mappings : IterableTypes
        Iterable of all dictionaries to be merged.
    collision_policy : optional[MergeCollisionPolicy]
        **Key collision policy** (i.e., strategy for merging keys shared by
        one or more mappings) to apply. Defaults to
        :attr:`MergeCollisionPolicy.RAISE_EXCEPTION`, raising an exception on
        the first key collision.
    is_values_copied : optional[bool]
        Either:

        * ``True``, if all values of the merged dictionary returned by this
          function are deeply (i.e., recursively) copied from their source
          dictionaries.
        * ``False``, if all values of the merged dictionary returned by this
          function are shallowly (i.e., non-recursively) copied from their
          source dictionaries.

        Defaults to ``False``.

    Returns
    ----------
    MappingType
        Dictionary merged from and of the same type as the passed dictionaries.

    Raises
    ----------
    BetseMappingException
        If the passed ``collision_policy`` is
        :attr:`MergeCollisionPolicy.RAISE_EXCEPTION` *and* any key of any
        key-value pair of any passed dictionary is also a key of any
        key-value pair of other such dictionary whose values differ.
    BetseSequenceException
        If less than two mappings are passed.

    See Also
    ----------
    http://treyhunner.com/2016/02/how-to-merge-dictionaries-in-python
        Blog post strongly inspiring this implementation. Thanks, Trey!
    '''

    # Avoid circular import dependencies.
    from betse.util.type.iterable import itertest, sequences
    from betse.util.type.iterable.mapping import maptest

    # Sequence of mappings converted from this iterable of mappings.
    mappings = sequences.to_sequence(mappings)

    # If less than two mappings were passed, raise an exception.
    sequences.die_if_length_less_than(sequence=mappings, length=2)
    # Else, at least two mappings were passed.

    # If any passed mapping is *NOT* a mapping, raise an exception.
    itertest.die_unless_items_instance_of(iterable=mappings, cls=MappingType)

    # Type of dictionary to be returned.
    dict_type = type(mappings[0])

    # Function accepting a single value and returning either a shallow or deep
    # copy of that value depending on which the caller requested.
    copy_value = deepcopy if is_values_copied else lambda value: value

    # If raising an exception on the first key collision, do so *BEFORE*
    # performing this merger.
    #
    # While there exist a countably infinite number of approaches to merging
    # non-colliding dictionaries in Python, this is the optimally efficient.
    # The is_maps_collide() function underlying the die_if_maps_collide()
    # function called below reduces to two symmetric set differences of
    # arbitrarily many iterables. Likewise, the merger performed below
    # transparently supports this key collision policy assuming no collisions.
    if collision_policy is MergeCollisionPolicy.RAISE_EXCEPTION:
        maptest.die_if_maps_collide(*mappings)
    # If giving higher precedence to dictionaries passed earlier, reverse the
    # order of the passed dictionaries. Why? Because the algorithm implemented
    # below implements the "PREFER_LAST" rather than "PREFER_FIRST" policy by
    # default. Since the former is merely the reverse of the latter, reversing
    # the order of these dictionaries repeatedly replaces the subsequent value
    # of each colliding key from the subsequent dictionary with the prior value
    # of that key in the prior dictionary, efficiently implementing the
    # "PREFER_FIRST" policy.
    elif collision_policy is MergeCollisionPolicy.PREFER_FIRST:
        mappings = reversed(mappings)

    # Dictionary merged from the passed dictionaries via a doubly-nested
    # dictionary comprehension (in the passed order). This repeatedly replaces
    # the prior value of each colliding key from the prior dictionary with the
    # subsequent value of that key in the subsequent dictionary, efficiently
    # implementing the "PREFER_LAST" policy.
    #
    # While there exist a countably infinite number of approaches to merging
    # dictionaries in Python, this approach is known to be the most efficient
    # for general-purpose merging of arbitrarily many dictionaries under Python
    # >= 3.4. See also Trey Hunter's exhaustive commentary (complete with
    # timings) at the above URL.
    dict_merged = {
        key: copy_value(value)
        for mapping in mappings
        for key, value in mapping.items()
    }

    # Return a dictionary of this type converted from this dictionary. If the
    # desired type is a "dict", this dictionary is returned as is; else, this
    # dictionary is converted into an instance of the desired type.
    return dict_merged if dict_type is dict else dict_type(dict_merged)
