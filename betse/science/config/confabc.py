#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Abstract base classes of all YAML-backed simulation configuration subclasses as
well as functionality pertaining to such classes.
'''

# ....................{ IMPORTS                            }....................
from abc import ABCMeta, abstractmethod
from betse.util.type.cls import classes
from betse.util.type.cls.descriptors import expr_alias, expr_enum_alias
from betse.util.type.obj import objects
from betse.util.type.types import (
    type_check, ClassType, EnumType, MappingType, SequenceTypes, TestableTypes)
from collections.abc import MutableSequence

# ....................{ SUPERCLASSES                       }....................
class SimConfABC(object, metaclass=ABCMeta):
    '''
    Abstract base class of all simulation configuration subclasses, each
    encapsulating a dictionary of related configuration settings (e.g.,
    representing one tissue profile) both loaded from and savable back to the
    current YAML-formatted simulation configuration file.

    Attributes
    ----------
    _conf : MappingType
        Low-level dictionary of related configuration settings both loaded from
        and savable back to the current YAML-formatted simulation configuration
        file.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(self, conf: MappingType) -> None:
        '''
        Associate this simulation configuration with the passed dictionary.

        Attributes
        ----------
        conf : MappingType
            Dictionary of related configuration settings both loaded from and
            savable back to the current YAML-formatted simulation configuration
            file.
        '''

        # Classify all passed parameters.
        self._conf = conf

    # ..................{ PROPERTIES ~ read-only             }..................
    # Read-only properties, preventing callers from resetting these attributes.

    @property
    def conf(self) -> MappingType:
        '''
        Dictionary of related configuration settings both loaded from and
        savable back to the current YAML-formatted simulation configuration
        file.
        '''

        return self._conf


class SimConfListableABC(SimConfABC):
    '''
    Abstract base class of all simulation configuration subclasses intended to
    be added as elements of lists of type :class:`SimConfList`.
    '''

    # ..................{ SUBCLASS                           }..................
    @abstractmethod
    def default(self) -> None:
        '''
        Initialize the dictionary associated with this simulation configuration
        to default values.

        This method is principally intended to be called by the
        :meth:`SimConfList.append_default` method, appending a new instance of
        this class (initialized by this method) to a list of these instances.
        '''

        pass

# ....................{ SUPERCLASSES ~ list                }....................
class SimConfList(MutableSequence):
    '''
    Simulation configuration list encapsulating a list of all dictionaries of
    related configuration settings (e.g., representing all tissue profiles) both
    loaded from and savable back to the current YAML-formatted simulation
    configuration file.

    Attributes
    ----------
    _confs_wrap : list
        High-level list of all instances of the :attr:`_conf_type` subclass
        encapsulating each dictionary in the low-level :attr:`_confs_yaml` list.
    _confs_yaml : SequenceTypes
        Low-level list of all dictionaries of related configuration settings
        both loaded from and savable back to the current YAML-formatted
        simulation configuration file.
    _conf_type : ClassType
        Subclass of the :class:`SimConfListableABC` abstract base class with which
        to instantiate each simulation configuration object encapsulating
        each dictionary in the :attr:`_confs_yaml` list.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(self, confs: SequenceTypes, conf_type: ClassType) -> None:
        '''
        Initialize this simulation configuration sublist.

        Attributes
        ----------
        confs : MappingType
            List of all dictionaries of related configuration settings both
            loaded from and savable back to the current YAML-formatted
            simulation configuration file.
        conf_type : ClassType
            Subclass of the :class:`SimConfListableABC` abstract base class with which
            to instantiate each simulation configuration object encapsulating
            each dictionary in the passed ``confs`` list. Specifically, for each
            such dictionary, a new object of this subclass is appended to the
            internal :attr:`_confs_wrap` list of such objects.
        '''

        # Raise an exception unless the passed type subclasses the desired API.
        classes.die_unless_subclass(subclass=conf_type, superclass=SimConfListableABC)

        # Classify all passed parameters.
        self._confs_yaml = confs
        self._conf_type = conf_type

        # Wrap each dictionary in this list with a new object of this type.
        self._confs_wrap = []
        for conf_yaml in self._confs_yaml:
            self._confs_wrap.append(self._conf_type(conf=conf_yaml))

    # ..................{ SUPERCLASS                         }..................
    # Abstract methods required by our superclass.

    def __len__(self):
        return len(self._confs_wrap)


    @type_check
    def __delitem__(self, index: int) -> None:
        del self._confs_yaml[index]
        del self._confs_wrap[index]


    @type_check
    def __getitem__(self, index: int) -> SimConfListableABC:
        return self._confs_wrap[index]


    @type_check
    def __setitem__(self, index: int, value: SimConfListableABC) -> None:
        '''
        Set the simulation configuration instance with the passed 0-based index
        to the passed such instance.
        '''

        # Raise an exception unless the passed object is an instance of the
        # desired API, which is specified at initialization time and hence
        # cannot be type checked above by a method annotation.
        objects.die_unless_instance(obj=value, cls=self._conf_type)

        # Set the high-level list item with this index to this object.
        self._confs_wrap[index] = value

        # Set the low-level list item with this index to this object's
        # underlying YAML-backed dictionary.
        self._confs_yaml[index] = value.conf


    @type_check
    def insert(self, index: int, value: SimConfListableABC) -> None:
        '''
        Insert the passed simulation configuration instance immediately *before*
        the simulation configuration instance with the passed 0-based index.
        '''

        # Raise an exception unless the passed object is an instance of the
        # desired API, which is specified at initialization time and hence
        # cannot be type checked above by a method annotation.
        objects.die_unless_instance(obj=value, cls=self._conf_type)

        # Insert this object *BEFORE* the high-level list item with this index.
        self._confs_wrap.insert(index, value)

        # Insert this object's underlying YAML-backed dictionary *BEFORE* the
        # low-level list item with this index.
        self._confs_yaml.insert(index, value.conf)

    # ..................{ APPENDERS                          }..................
    def append_default(self) -> SimConfListableABC:
        '''
        Append a new simulation configuration instance initialized to default
        values immediately *after* the last such instance.

        Returns
        ----------
        SimConfListableABC
            Instance of the :attr:`_conf_type` subclass appended to the
            high-level :attr:`_confs_wrap` list, encapsulating the new
            dictionary appended to the low-level :attr:`_confs_yaml` list.
        '''

        # YAML-backed dictionary underlying the new instance to be appended.
        conf_yaml = {}

        # Instance wrapping this dictionary.
        conf_wrap = self._conf_type(conf=conf_yaml)

        # Initialize this instance to default values.
        conf_wrap.default()

        # Return this instance.
        return conf_wrap

# ....................{ DESCRIPTORS                        }....................
@type_check
def conf_alias(keys: str, cls: TestableTypes = None) -> object:
    '''
    Expression alias **data descriptor** (i.e., object satisfying the data
    descriptor protocol) specific to simulation configurations, dynamically
    aliasing a target variable bound to instances of the class instantiating
    this descriptor to a source Python expression performing one or more key
    lookups into the dictionary loaded from a YAML-formatted simulation
    configuration file.

    Parameters
    ----------
    keys : str
        Python expression evaluating to the value of an arbitrarily nested key
        of the dictionary loaded from the current simulation configuration,
        typically consisting of one or more ``[``- and ``]``-delimited key
        lookups into this same dictionary (e.g.,
        ``['variable settings']['noise']['dynamic noise']``).
    cls: optional[ClassType]
        Either the expected type or tuple of the expected types of this variable
        *or* a callable validating the values to which this variable is to be
        set. Defaults to ``None``, in which case this variable is permissively
        settable to *any* values.

    Returns
    ----------
    object
        Expression alias data descriptor as detailed above.

    See Also
    ----------
    :func:`expr_alias`
        Further details.
    '''

    return expr_alias(expr='self._conf' + keys, cls=cls)


@type_check
def conf_enum_alias(keys: str, enum_type: EnumType) -> object:
    '''
    Enumeration-specific expression alias **data descriptor** (i.e., object
    satisfying the data descriptor protocol) specific to simulation
    configurations, dynamically aliasing a target variable of the passed
    enumeration type bound to instances of the class instantiating this
    descriptor to an arbitrarily complex source Python expression performing one
    or more key lookups into the dictionary loaded from a YAML-formatted
    simulation configuration file.

    Parameters
    ----------
    keys : str
        Oner or more ``[``- and ``]``-delimited key lookups. See the
        :func:`conf_alias` function for further details.
    enum_type: EnumType
        Enumeration that the value of this variable *must* be a member of.
        Setting this variable to a value *not* a member of this enumeration will
        raise an exception.

    Returns
    ----------
    object
        Enumeration-specific expression alias data descriptor as detailed above.

    See Also
    ----------
    :func:`expr_enum_alias`
        Further details.
    '''

    return expr_enum_alias(expr='self._conf' + keys, enum_type=enum_type)
