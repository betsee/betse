#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level **mapping utilities** (i.e., functions operating on dictionary-like
types and instances).
'''

# ....................{ IMPORTS                            }....................
import pprint
from betse.exceptions import BetseMappingException
from betse.util.type.types import (
    type_check, MappingType, HashableType,)

# ....................{ EXCEPTIONS                         }....................
@type_check
def die_unless_values_unique(mapping: MappingType) -> None:
    '''
    Raise an exception unless *all* values of the passed dictionary are
    **unique** (i.e., no two values of two distinct key-value pairs are equal).

    Parameters
    ----------
    mapping : MappingType
        Dictionary to be inspected.

    Raises
    ----------
    BetseMappingException
        If at least one value of this dictionary is a duplicate.
    '''

    # If one or more values of this dictionary are duplicates...
    if not is_values_unique(mapping):
        # Avoid circular import dependencies.
        from betse.util.type import iterables

        # Set of all duplicate values in this dictionary.
        values_duplicate = iterables.get_items_duplicate(mapping.values())

        # Raise an exception embedding this set.
        raise BetseMappingException(
            'Dictionary values {} duplicate.'.format(values_duplicate))

# ....................{ TESTERS                            }....................
@type_check
def is_keys(mapping: MappingType, *keys: HashableType) -> bool:
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

    # Yes, this is ridiculously awesome.
    return set(keys).issubset(mapping)


@type_check
def is_values_unique(mapping: MappingType) -> bool:
    '''
    ``True`` only if *all* values of the passed dictionary are **unique** (i.e.,
    no two values of two distinct key-value pairs are equal).

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
    from betse.util.type import iterables

    # For sanity, defer to an existing low-level tester.
    return iterables.is_items_unique(mapping.values())

# ....................{ FORMATTERS                         }....................
@type_check
def format(mapping: MappingType) -> str:
    '''
    Convert the passed dictionary into a human-readable string.
    '''

    return pprint.pformat(mapping)

# ....................{ MERGERS                            }....................
@type_check
def merge(*dicts: MappingType) -> MappingType:
    '''
    Dictionary of all key-value pairs merged together from all passed
    dictionaries (in the passed order).

    Parameters
    ----------
    mappings : Tuple[MappingType]
        Tuple of all dictionaries to be merged.

    Returns
    ----------
    MappingType
        Dictionary merged from and of the same type as the passed dictionaries.
        For efficiency, this dictionary is only a shallow rather than deep copy
        of these dictionaries. Note lastly that the class of the passed
        dictionary *must* define an ``__init__()`` method accepting a dictionary
        comprehension.

    See Also
    ----------
    :meth:`dict.update`
        Standard method merging two dictionaries, which should typically be
        called instead of this slower function in this specific use case.
    '''

    # Type of dictionary to be returned.
    dict_type = type(dicts[0])

    # Dictionary merged from the passed dictionaries via a doubly-nested
    # dictionary comprehension. While there exist a countably infinite number of
    # approaches to merging dictionaries in Python, this approach is known to be
    # the most efficient for general-purpose merging of arbitrarily many
    # dictionaries under Python >= 3.4. See also Trey Hunter's exhaustive
    # commentary replete with timings at:
    #     http://treyhunner.com/2016/02/how-to-merge-dictionaries-in-python
    dict_merged = {
        key: value
        for dict_cur in dicts
        for key, value in dict_cur.items()
    }

    # Return a dictionary of this type converted from this dictionary. If the
    # desired type is a "dict", this dictionary is returned as is; else, this
    # dictionary is converted into an instance of the desired type.
    return dict_merged if dict_type is dict else dict_type(dict_merged)
