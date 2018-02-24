#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level **mapping classes** (i.e., classes implementing dictionary-like
functionality, typically by subclassing the builtin :class:`dict` container type
or an analogue thereof).
'''

# ....................{ IMPORTS                            }....................
from betse.exceptions import (
    BetseMappingException, BetseMethodUnimplementedException)
from betse.util.type import types
from betse.util.type.types import (
    type_check,
    CallableTypes,
    HashableType,
    IteratorType,
    MappingType,
    MappingMutableType,
)
from collections import OrderedDict

# ....................{ GLOBALS                            }....................
_DEFAULT_DICT_ID = 0
'''
Unique arbitrary identifier with which to uniquify the class name of the next
:func:`DefaultDict`-derived type.
'''

# ....................{ CLASSES ~ default                  }....................
#FIXME: Donate back to StackOverflow. The standard "defaultdict" class is
#sufficiently useless that numerous users would probably find this useful.
@type_check
def DefaultDict(
    #FIXME: For disambiguity, rename to "missing_key_maker".
    missing_key_value: CallableTypes,
    initial_mapping: MappingType = None,
) -> MappingType:
    '''
    Sane **default dictionary** (i.e., dictionary implicitly setting an unset
    key to the value returned by a caller-defined callable passed both this
    dictionary and that key).

    Motivation
    ----------
    The standard :class:`collections.defaultdict` class is sadly insane.
    Happily, this custom class is not. Specifically, that class:

    * Requires that the caller-defined callable accept *no* arguments.
    * Requires that the initial key-value pairs to seed this dictionary with be
      awkwardly defined as keyword arguments rather than key-value pairs.

    Parameters
    ----------
    missing_key_value : CallableTypes
        Callable (e.g., function, lambda, method) called to generate the default
        value for a "missing" (i.e., undefined) key on the first attempt to
        access that key, passed first this dictionary and then this key and
        returning this value. This callable should have a signature resembling:
        ``def missing_key_value(self: DefaultDict, missing_key: object) ->
        object``. Equivalently, this callable should have the exact same
        signature as that of the optional :meth:`dict.__missing__` method.
    initial_mapping : optional[MappingType]
        Non-default dictionary providing the initial key-value pairs of this
        default dictionary. Defaults to ``None``, in which case this
        default dictionary is initially empty.

    Returns
    ----------
    MappingType
        Default dictionary initialized as described above.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.cls import classes

    # Dynamically generated default dictionary class specific to this callable.
    default_dict_class = classes.define_class(
        class_name=_get_default_dict_class_name(),
        class_attr_name_to_value={'__missing__': missing_key_value,},
        base_classes=(dict,),
    )

    # Instantiate the first and only instance of this class.
    default_dict = default_dict_class()

    # If an initial dictionary was passed, add all of its key-value pairs.
    if initial_mapping is not None:
        default_dict.update(initial_mapping)

    # Return this instance.
    return default_dict

# ....................{ CLASSES ~ dynamic                  }....................
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


class DynamicValueDict(MappingMutableType):
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
    _key_to_dynamic_value : MappingType
        Dictionary such that each:
        * Key is a key of the outer dictionary.
        * Value is the `DynamicValue` instance permitting that key's underlying
          variable to be dynamically get and set.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(self, key_to_dynamic_value: MappingType) -> None:
        '''
        Initialize this dictionary.

        Parameters
        ----------
        key_to_dynamic_value : MappingType
            Dictionary such that each:
            * Key is a key of the outer dictionary.
            * Value is the `DynamicValue` instance permitting that key's
              underlying variable to be dynamically get and set.
        '''

        #FIXME: Premature optimization. Reduce this to the following, please:
        #    from betse.util.type import iterables
        #    iterables.die_unless_items_instance_of(
        #        iterable=key_to_dynamic_value.values(), cls=DynamicValue)

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

# ....................{ CLASSES ~ ordered                  }....................
class OrderedArgsDict(OrderedDict):
    '''
    Ordered dictionary initialized by a sequence of key-value pairs.

    Motivation
    ----------
    This :class:`OrderedDict` subclass simplifies the definition of ordered
    dictionaries. As is common practice throughout the stdlib, the standard
    :class:`OrderedDict` class is defined by a sequence of 2-sequences
    ``(key, value)`` (e.g., tuple of 2-tuples) where each ``key`` and ``value``
    comprises one key-value pair of this dictionary.

    This :class:`OrderedDict` subclass is instead defined by a flat sequence of
    keys and values rather than a nested sequence of 2-sequences. In all other
    respects, this subclass is identical to the :class:`OrderedDict` class.

    Examples
    ----------
    >>> from betse.util.type.mapping.mapcls import OrderedArgsDict
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
        from betse.util.type.numeric import ints

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

# ....................{ CLASSES ~ reversible               }....................
#FIXME: Define an analogous "OneToOneDict" subclass of the "dict" superclass,
#mandating that no values ambiguously map to by two or more keys of the original
#dictionary. In theory, the implementation should be even more concise.
#FIXME: Unit test us up.
class ReversibleDict(dict):
    '''
    Bidirectional dictionary.

    This dictionary subclass provides a public :attr:``reverse`` instance
    variable, whose value is a dictionary providing both safe and efficient
    **reverse lookup** (i.e., lookup by values rather than keys) of the
    key-value pairs of the original dictionary.

    This dictionary explicitly supports reverse lookup of values ambiguously
    mapped to by two or more keys of the original dictionary. To disambiguate
    between these mappings, each key of the :attr:``reverse`` dictionary
    (signifying a value of the original dictionary) maps to a sequence
    containing all keys of the original dictionary mapping to that value.

    Caveats
    ----------
    **The :attr:``reverse`` dictionary is not safely modifiable.**

    Attempting to do so breaks the contractual semantics of this class in an
    unsupported manner. For safety, this dictionary is silently synchronized
    with *all* modifications to the original dictionary -- including additions
    of new keys, modifications of the values of existing keys, and deletions of
    existing keys. However, the reverse is *not* the case; the original
    dictionary is *not* silently synchronized with any modifications to the
    :attr:``reverse`` dictionary. (Doing so is technically feasible but beyond
    the scope of current use cases.)

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

# ....................{ PRIVATE                            }....................
def _get_default_dict_class_name() -> str:
    '''
    Name of the class of the next data descriptor created and returned by the
    next call to the :func:`DefaultDict` function, guaranteed to be unique
    across all such classes.

    Since the human-readability of this name is neither required nor desired,
    this is a non-human-readable name efficiently constructed from a private
    prefix and a unique arbitrary identifier.
    '''

    # Global variable modified below.
    global _DEFAULT_DICT_ID

    # Increment this identifier to preserve uniqueness.
    _DEFAULT_DICT_ID += 1

    # Return a unique classname suffixed by this identifier.
    return '__DefaultDict' + str(_DEFAULT_DICT_ID)
