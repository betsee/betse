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
from copy import deepcopy

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

    # Avoid circular import dependencies.
    from betse.util.type import iterables
    from betse.util.type.text import strs

    # If one or more values of this dictionary are duplicates...
    if not is_values_unique(mapping):
        # Set of all duplicate values in this dictionary.
        values_duplicate = iterables.get_items_duplicate(mapping.values())

        # Raise an exception embedding this set.
        raise BetseMappingException(
            'Dictionary values {} duplicate.'.format(
                strs.join_as_conjunction_double_quoted(*values_duplicate)))

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

# ....................{ COPIERS                            }....................
@type_check
def copy(mapping: MappingType) -> MappingType:
    '''
    Dictionary of all key-value pairs deeply (i.e., recursively) duplicated from
    the passed dictionary.

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

# ....................{ MERGERS                            }....................
@type_check
def merge(*dicts: MappingType) -> MappingType:
    '''
    Dictionary of all key-value pairs deeply (i.e., recursively) merged together
    from all passed dictionaries (in the passed order).

    **Order is significant.** Dictionaries passed later take precedence over
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
        # For safety, deeply copy rather than reuse this value.
        key: deepcopy(value)
        for dict_cur in dicts
        for key, value in dict_cur.items()
    }

    # Return a dictionary of this type converted from this dictionary. If the
    # desired type is a "dict", this dictionary is returned as is; else, this
    # dictionary is converted into an instance of the desired type.
    return dict_merged if dict_type is dict else dict_type(dict_merged)
