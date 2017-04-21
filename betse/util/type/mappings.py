#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level dictionary facilities.
'''

# ....................{ IMPORTS                            }....................
import pprint
from betse.exceptions import (
    BetseMappingException, BetseMethodUnimplementedException)
from betse.util.type import types
from betse.util.type.types import (
    type_check,
    CallableTypes,
    HashableType,
    IteratorType,
    MappingType,
    MutableMappingType,
)
from collections import OrderedDict

# ....................{ CLASSES ~ dict : dynamic           }....................
class DynamicValue(object):
    '''
    **Dynamic value** (i.e., qbject encapsulating a single variable gettable and
    settable via callables predefined at object initialization time).

    This object typically encapsulates a variable repeatedly reassigned to.
    Since the value of such a variable is inconstant, this variable is *not*
    reliably passable or referrable to without introducing desynchronization
    issues (e.g., the ``sim.cc_cells[sim.iNa]`` Numpy array of all cell-specific
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

    # ..................{ SLOTS                              }..................
    # Tuple of the names of all instance attributes permitted in instances of
    # this class. This slightly improves the time efficiency of attribute access
    # (by anywhere from 5% to 10%) and dramatically improves the space
    # efficiency of object storage (by several orders of magnitude).
    __slots__ = ('get_value', 'set_value',)

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(
        self, get_value: CallableTypes, set_value: CallableTypes) -> None:
        '''
        Initialize this dynamic value.

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

    # ..................{ INITIALIZERS                       }..................
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

# ....................{ CLASSES ~ dict : ordered           }....................
class OrderedArgsDict(OrderedDict):
    '''
    Ordered dictionary initialized by a sequence of key-value pairs.

    The canonical :class:`OrderedDict` class implements ordered dictionaries in
    Python. As is common practice throughout the stdlib, this class is
    initialized with a sequence of 2-sequences (e.g., tuple of 2-tuples) whose
    first elements comprise the keys and whose second elements comprise the
    values of all key-value pairs of this dictionary.

    This :class:`OrderedDict` subclass simplifies this initializing burden.
    Specifically, this class is initialized with a flat sequence of keys and
    values rather than a nested sequence of 2-sequences. In all other respects,
    this class is identical to the :class:`OrderedDict` class.

    Examples
    ----------
    >>> from betse.util.type.mappings import OrderedArgsDict
    >>> from collections import OrderedDict
    >>> tuatha_de_danann = OrderedArgsDict(
    ...     'Nuada', 'Nodens',
    ...     'Lugh', 'Lugus',
    ...     'Brigit', 'Brigantia',
    ... )
    >>> tuath_de = OrderedDict((
    ...     ('Nuada', 'Nodens'),
    ...     ('Lugh', 'Lugus'),
    ...     ('Brigit', 'Brigantia'),
    ... ))
    >>> for irish_name, celtic_name in tuatha_de_danann.items():
    ...     print('From Irish {} to Celtic {}.'.format(irish_name, celtic_name))
    From Irish Nuada to Celtic Nodens.
    From Irish Lugh to Celtic Lugus.
    From Irish Brigit to Celtic Brigantia.
    >>> tuatha_de_danann == tuath_de
    True
    '''


    def __init__(self, *key_value_pairs) -> None:
        '''
        Initialize this ordered dictionary with the passed tuple of sequential
        keys and values of even length.

        Attributes
        ----------
        key_value_pairs : tuple
            Tuple of sequential keys and values to initialize this ordered
            dictionary with. Each element of this tuple with:
            * Even index (e.g., the first and third elements) defines a new key
              of this dictionary, beginning a new key-value pair.
            * Odd index (e.g., the second and fourth elements) defines the value
              for the key defined by the preceding element, finalizing this
              existing key-value pair.

        Raises
        ----------
        :exc:`betse.exceptions.BetseMappingException`
             If this tuple is _not_ of even length.
        '''

        # Avoid circular import dependencies.
        from betse.util.type import ints

        # If the passed tuple of key-value pairs is odd and hence omitted the
        # value for the final key, raise an exception.
        if ints.is_odd(len(key_value_pairs)):
            raise BetseMappingException(
                'Expected even number of key-value parameters, '
                'but received {} such parameters.'.format(len(key_value_pairs)))

        # Zip object yielding 2-tuple key-value pairs of the format required by
        # the superclass, converted from this flat sequence of keys and values.
        key_value_pairs_nested = zip(
            key_value_pairs[0::2],
            key_value_pairs[1::2])

        # Initialize the superclass with this zip.
        super().__init__(key_value_pairs_nested)

# ....................{ CLASSES ~ dict : reversible        }....................
#FIXME: Unit test this.
class ReversibleDict(dict):
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

# ....................{ TESTERS                            }....................
@type_check
def is_keys(mapping: MappingType, *keys: HashableType) -> bool:
    '''
    `True` only if the passed dictionary contains _all_ passed keys.

    Parameters
    ----------
    mapping : MappingType
        Dictionary to be tested.
    keys : tuple[HashableType]
        Tuple of all keys to be tested for.

    Returns
    ----------
    bool
        `True` only if this dictionary contains _all_ passed keys.
    '''

    # Yes, this is ridiculously awesome.
    return set(keys).issubset(mapping)

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
    dictionaries (_in the passed order_).

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
        dictionary _must_ define an `__init__()` method accepting a dictionary
        comprehension.
    '''

    # Type of dictionary to be returned.
    dict_type = type(dicts[0])

    # Dictionary merged from the passed dictionaries via a doubly-nested
    # dictionary comprehension. While there exist a countably infinite number of
    # approaches to merging dictionaries in Python, this approach is known to be
    # the most efficient for general-purpose merging of arbitrarily many
    # dictionaries under Python >= 3.4. See also Trey Hunter's exhaustive
    # commentary replete with timings at:
    #     http://treyhunner.com/2016/02/how-to-merge-dictionaries-in-python/
    dict_merged = {
        key: value
        for dict_cur in dicts
        for key, value in dict_cur.items()
    }

    # Return a dictionary of this type converted from this dictionary. If the
    # desired type is a "dict", this dictionary is returned as is; else, this
    # dictionary is converted into an instance of the desired type.
    return dict_merged if dict_type is dict else dict_type(dict_merged)
