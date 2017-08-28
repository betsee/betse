#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Abstract base classes of all YAML-backed simulation configuration subclasses
encapsulating on-disk lists and list items.
'''

# ....................{ IMPORTS                            }....................
from abc import abstractmethod
from betse.lib.yaml.yamlabc import YamlABC
from betse.lib.yaml.yamlalias import yaml_alias
from betse.science.simulate.pipe.piperun import SimPipeRunnerConfMixin
from betse.util.type.cls import classes
from betse.util.type.obj import objects
from betse.util.type.types import type_check, SequenceOrNoneTypes, ClassType
from collections import MutableSequence

# ....................{ SUPERCLASSES ~ list item           }....................
class YamlListItemABC(SimPipeRunnerConfMixin, YamlABC):
    '''
    Abstract base class of all simulation list item subconfigurations, each
    backed by a YAML list item and intended to be added to a
    :class:`YamlList` container.

    Design
    ----------
    This class subclasses the :class:`SimPipeRunnerConfMixin` mixin, allowing
    all instances of:

    * This class to be used as **simulation pipeline runner arguments** (i.e.,
      simplistic objects encapsulating all input parameters passed to runner
      methods in :class:`SimPipeABC` pipelines).
    * The :class:`YamlList` class to be used as sequences of these arguments
      and hence returned from the abstract
      :class:`SimPipeABC._runners_conf_enabled` property.
    '''

    # ..................{ MAKERS                             }..................
    @classmethod
    @abstractmethod
    def make_default(cls) -> 'betse.science.config.confabc.YamlListItemABC':
        '''
        Create and return an instance of this subclass encapsulating a new
        dictionary containing default configuration settings.

        This method is principally intended to be called by the
        :meth:`YamlList.append_default` method, appending this instance to an
        existing list of such instances.
        '''

        pass


class YamlListItemTypedABC(YamlListItemABC):
    '''
    Abstract base class of all simulation typed list item subconfigurations,
    each backed by a YAML list item whose dictionary keys define the type of
    this item and intended to be added to a :class:`YamlList` container.

    Attributes
    ----------
    is_enabled : bool
        ``True`` only if this list item is enabled.
    name : str
        Lowercase alphanumeric string uniquely identifying the type of this
        list item (e.g., ``voltage_membrane``, signifying a transmembrane
        voltage list item). See each ``type`` key of the corresponding list in
        the default simulation configuration file for further commentary.
    '''

    # ..................{ ALIASES                            }..................
    is_enabled = yaml_alias("['enabled']", bool)
    name       = yaml_alias("['type']", str)

# ....................{ SUPERCLASSES ~ list                }....................
class YamlList(MutableSequence):
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
        Subclass of the :class:`YamlListItemABC` abstract base class with
        which to instantiate each simulation configuration object encapsulating
        each dictionary in the :attr:`_confs_yaml` list.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(
        self, confs: SequenceOrNoneTypes, conf_type: ClassType) -> None:
        '''
        Initialize this simulation configuration sublist.

        Attributes
        ----------
        confs : MappingType
            List of all dictionaries of related configuration settings both
            loaded from and savable back to the current YAML-formatted
            simulation configuration file *or* ``None`` if the key defining this
            list in this file has no corresponding value, in which case this
            list defaults to the empty list.
        conf_type : ClassType
            Subclass of the :class:`YamlListItemABC` abstract base class
            with which to instantiate each simulation configuration object
            encapsulating each dictionary in the passed ``confs`` list.
            Specifically, for each such dictionary, a new object of this type
            is appended to the internal :attr:`_confs_wrap` list of these
            objects.
        '''

        # Raise an exception unless the passed type implements the listable API.
        classes.die_unless_subclass(
            subclass=conf_type, superclass=YamlListItemABC)

        # If this list is unspecified, default this list to the empty list.
        if confs is None:
            confs = []

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
    def __getitem__(self, index: int) -> YamlListItemABC:
        return self._confs_wrap[index]


    @type_check
    def __setitem__(self, index: int, value: YamlListItemABC) -> None:
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
    def insert(self, index: int, value: YamlListItemABC) -> None:
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
    def append_default(self) -> YamlListItemABC:
        '''
        Append a new simulation configuration list item initialized to default
        values to this list and return this list item.

        Returns
        ----------
        YamlListItemABC
            Instance of the :attr:`_conf_type` subclass appended to the
            high-level :attr:`_confs_wrap` list, encapsulating the new
            dictionary appended to the low-level :attr:`_confs_yaml` list.
        '''

        # Default simulation configuration list item.
        conf_wrap = self._conf_type.make_default()

        # Append this list item to this list.
        self.append(conf_wrap)

        # Return this list item.
        return conf_wrap
