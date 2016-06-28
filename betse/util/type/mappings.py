#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level dictionary facilities.
'''

# ....................{ IMPORTS                            }....................
import pprint
from betse.util.type.types import type_check
from collections.abc import Mapping

# ....................{ CLASSES                            }....................
class bidict(dict):
    '''
    Bidirectional dictionary.

    This dictionary subclass provides a public `reverse` attribute, whose value
    is a dictionary providing both safe and efficient **reverse lookup** (i.e.,
    lookup by values rather than keys) of the key-value pairs of the original
    dictionary.

    The `reverse` dictionary explicitly supports reverse lookup of values
    ambiguously mapped to by two or more keys of the original dictionary. To
    disambiguate between these mappings, each key of the `reverse` dictionary
    (signifying a value of the original dictionary) maps to a list
    containing all keys of the original dictionary mapping to that value.

    Caveats
    ----------
    **The `reverse` dictionary is _not_ safely modifiable.** Attempting to do so
    breaks the contractual semantics of this class in an unsupported manner. For
    safety, the `reverse` dictionary is silently synchronized with _all_
    modifications to the original dictionary -- including additions of new keys,
    modifications of the values of existing keys, and deletions of existing
    keys. However, the reverse is _not_ the case; the original dictionary is
    _not_ silently synchronized with any modifications to the `reverse`
    dictionary. (Doing so is technically feasible but beyond the scope of this
    class' use, currently.)

    Attributes
    ----------
    reverse : dict
        Dictionary such that each:
        * Key is a value of the original dictionary.
        * Value is a list of all keys of the original dictionary mapping to that
          value of the original dictionary.

    See Also
    ----------
    https://stackoverflow.com/a/21894086/2809027
        Stackoverflow answer strongly inspiring this class. Thanks alot, Basj!
    '''


    def __init__(self, *args, **kwargs):
        # Initialize the original dictionary.
        super().__init__(*args, **kwargs)

        # Initialize the reversed dictionary. For each key-value pair with which
        # the original dictionary is initialized, map this value to a list of
        # all corresponding keys for reverse lookup.
        self.reverse = {}
        for key, value in self.items():
            self.reverse.setdefault(value, []).append(key)


    def __setitem__(self, key, value):
        # Map this key to this value in the original dictionary.
        super().__setitem__(key, value)

        # Map this value to this key in the reversed dictionary.
        self.reverse.setdefault(value, []).append(key)


    def __delitem__(self, key):
        # Value to be deleted from the reversed dictionary.
        value = self[key]

        # Remove this key-value pair from the original dictionary.
        super().__delitem__(key)

        # Remove this value-key pair from the reversed dictionary.
        self.reverse[value].remove(key)

        # If this key is the last key mapping to this value in the original
        # dictionary, remove this value entirely from the reversed dictionary.
        # While arguably optional, doing so reduces space costs with no
        # concomitant increase in time costs.
        if not self.reverse[value]:
            del self.reverse[value]

# ....................{ TESTERS                            }....................
@type_check
def is_keys(mapping: Mapping, *key_names) -> str:
    '''
    `True` only if the passed dictionary contains _all_ passed keys.

    Parameters
    ----------
    mapping : Mapping
        Dictionary to be tested.
    key_names : tuple
        Tuple of all keys to be tested for.
    '''

    # Yes, this is ridiculously awesome.
    return set(key_names).issubset(mapping)

# ....................{ FORMATTERS                         }....................
def format(mapping: Mapping) -> str:
    '''
    Convert the passed dictionary into a human-readable string.
    '''

    return pprint.pformat(mapping)
