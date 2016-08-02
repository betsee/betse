#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level dictionary facilities.
'''

# ....................{ IMPORTS                            }....................
import pprint
from betse.exceptions import BetseMethodUnimplementedException
from betse.util.type import types
from betse.util.type.types import (
    type_check,
    CallableTypes,
    HashableType,
    IteratorType,
    MappingType,
    MutableMappingType,
)

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
        Stackoverflow answer inspiring this class. Cheers, Basj!
    '''


    def __init__(self, *args, **kwargs) -> None:
        # Initialize the original dictionary.
        super().__init__(*args, **kwargs)

        # Initialize the reversed dictionary. For each key-value pair with which
        # the original dictionary is initialized, map this value to a list of
        # all corresponding keys for reverse lookup.
        self.reverse = {}
        for key, value in self.items():
            self.reverse.setdefault(value, []).append(key)


    def __setitem__(self, key: HashableType, value: object) -> None:
        # Map this key to this value in the original dictionary.
        super().__setitem__(key, value)

        # Map this value to this key in the reversed dictionary.
        self.reverse.setdefault(value, []).append(key)


    def __delitem__(self, key: HashableType) -> None:
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

# ....................{ CLASSES ~ dynamic                  }....................
class DynamicValue(object):
    '''
    Object encapsulating a single variable gettable and settable via callables
    predefined at object initialization time.

    This object typically encapsulates a variable repeatedly reassigned to.
    Since the value of such a variable is inconstant, this variable is _not_
    reliably passable or referrable to without introducing desynchronization
    issues (e.g., the `sim.cc_cells[sim.iNa]` Numpy array of all cell-specific
    sodium ion concentrations, unreliably reassigned to each simulation time
    step). This object wraps this variable's access and modification, permitting
    this variable to effectively be reliably passed and referred to without
    introducing such issues.

    Attributes
    ----------
    get_value : CallableTypes
        Callable (e.g., function, lambda, method) dynamically getting this
        variable's value.
    set_value : CallableTypes
        Callable (e.g., function, lambda, method) dynamically setting this
        variable's value.
    '''

    @type_check
    def __init__(
        self, get_value: CallableTypes, set_value: CallableTypes) -> None:
        '''
        Initialize this object.

        Parameters
        ----------
        get_value : CallableTypes
            Callable (e.g., function, lambda, method) dynamically getting this
            variable's value.
        set_value : CallableTypes
            Callable (e.g., function, lambda, method) dynamically setting this
            variable's value.
        '''

        self.get_value = get_value
        self.set_value = set_value


class DynamicValueDict(MutableMappingType):
    '''
    Dictionary whose values are all **dynamic** (i.e., instances of the
    `DynamicValue` class, whose actual underlying values are gettable and
    settable via callables specific to those instances).

    Caveats
    ----------
    **All dictionary keys are predefined at dictionary initialization time.**
    Unlike standard dictionaries, this dictionary implementation does _not_
    permit additional key-value pairs to be added to or existing key-value pairs
    to be removed from this dictionary after initialization. While the value to
    which any key maps is arbitrarily modifiable, no key itself is modifiable.

    Attributes
    ----------
    _key_to_dynamic_value : dict
        Dictionary such that each:
        * Key is a key of the outer dictionary.
        * Value is the `DynamicValue` instance permitting that key's underlying
          variable to be dynamically get and set.
    '''

    # ..................{ LIFECYCLE                          }..................
    @type_check
    def __init__(self, key_to_dynamic_value: dict) -> None:
        '''
        Initialize this dictionary.

        Parameters
        ----------
        key_to_dynamic_value : dict
            Dictionary such that each:
            * Key is a key of the outer dictionary.
            * Value is the `DynamicValue` instance permitting that key's
              underlying variable to be dynamically get and set.
        '''

        # If optimization is disabled, validate all values of the passed
        # dictionary to be instances of the "DynamicValue" class.
        if __debug__:
            for dynamic_value in key_to_dynamic_value.values():
                assert isinstance(dynamic_value, DynamicValue), (
                    '"{}" not a dynamic value.'.format(
                        types.trim(dynamic_value)))

        # Classify the passed parameters.
        self._key_to_dynamic_value = key_to_dynamic_value


    # ..................{ MAGIC ~ dict                       }..................
    def __getitem__(self, key: HashableType) -> object:
        '''
        Get the underlying value of the variable associated with the passed key.

        Parameters
        ----------
        key : HashableType
            Key to get the value of.

        Returns
        ----------
        object
            Underlying value of the variable associated with this key.

        Raises
        ----------
        KeyError
            If this key is _not_ a key with which this dictionary was
            initialized.
        '''

        return self._key_to_dynamic_value[key].get_value()


    def __setitem__(self, key: HashableType, value: object) -> None:
        '''
        Set the underlying value of the variable associated with the passed key.

        Parameters
        ----------
        key : HashableType
            Key to set the value of.
        value : object
            Underlying value of the variable associated with this key to be set.

        Raises
        ----------
        KeyError
            If this key is _not_ a key with which this dictionary was
            initialized.
        '''

        # Attempt to get this key *BEFORE* setting this key, implicitly raising
        # a "KeyError" exception if this key was not defined at initialization.
        # self._key_to_dynamic_value[key]

        # Set this key.
        self._key_to_dynamic_value[key].set_value(value)


    def __delitem__(self, key: HashableType) -> None:
        '''
        Unconditionally raise an exception.

        Since the size of this dictionary is invariant, deleting key-value pairs
        from this dictionary is strictly prohibited.

        Raises
        ----------
        BetseMethodUnimplementedException
            If this key is _not_ a key with which this dictionary was
            initialized.
        '''

        raise BetseMethodUnimplementedException()

    # ..................{ MAGIC ~ container                  }..................
    def __len__(self) -> int:
        '''
        Get the size of this dictionary, guaranteed to be equal to the size of
        the dictionary with which this dictionary was initialized.

        Returns
        ----------
        int
            This size.
        '''

        return len(self._key_to_dynamic_value)

    # ..................{ MAGIC ~ iterable                   }..................
    def __iter__(self) -> IteratorType:
        '''
        Get an iterator over all keys with which this dictionary was initalized.

        Returns
        ----------
        IteratorType
            This iterator.
        '''

        return iter(self._key_to_dynamic_value)

# ....................{ TESTERS                            }....................
@type_check
def is_keys(mapping: MappingType, *key_names) -> str:
    '''
    `True` only if the passed dictionary contains _all_ passed keys.

    Parameters
    ----------
    mapping : MappingType
        Dictionary to be tested.
    key_names : tuple
        Tuple of all keys to be tested for.
    '''

    # Yes, this is ridiculously awesome.
    return set(key_names).issubset(mapping)

# ....................{ FORMATTERS                         }....................
def format(mapping: MappingType) -> str:
    '''
    Convert the passed dictionary into a human-readable string.
    '''

    return pprint.pformat(mapping)
