#!/usr/bin/env python3
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Abstract base classes of all YAML-backed simulation configuration subclasses
encapsulating on-disk lists and list items.
'''

# ....................{ IMPORTS                            }....................
from abc import abstractmethod
from betse.lib.yaml.yamlalias import yaml_alias
from betse.lib.yaml.abc.yamlabc import YamlABC
from betse.util.type.cls import classes
from betse.util.type.obj import objects
from betse.util.type.types import (
    type_check, ClassType, SequenceTypes)
from collections import MutableSequence

# ....................{ SUPERCLASSES ~ list item           }....................
class YamlListItemABC(YamlABC):
    '''
    Abstract base class of all YAML-backed list item subclasses, each instance
    of which is intended to be added to a parent :class:`YamlList` container.

    Each such instance may technically encapsulate any valid YAML type (e.g.,
    :class:`int`, :class:`str`) but typically encapsulates a YAML dictionary of
    related key-value pairs (e.g., animation settings, tissue profile).
    '''

    # ..................{ MAKERS                             }..................
    @classmethod
    def make_list(cls, *args, **kwargs) -> (
        'betse.lib.yaml.abc.yamllistabc.YamlList'):
        '''
        Create and return a low-level YAML-backed list containing only items of
        this specific subclass type, both loaded from and savable back into a
        YAML-formatted file.

        Parameters
        ----------
        All parameters are passed as is to the :meth:`YamlList.__init__` method.
        '''

        return YamlList(*args, item_type=cls, **kwargs)

    # ..................{ MAKERS ~ abstract                  }..................
    # Subclasses are required to implement the following abstract class methods.

    @classmethod
    @abstractmethod
    def make_default(
        cls,
        yaml_list: 'betse.lib.yaml.abc.yamllistabc.YamlList',
    ) -> 'betse.lib.yaml.abc.yamllistabc.YamlListItemABC':
        '''
        Create and return an instance of this subclass encapsulating a new
        dictionary containing default configuration settings.

        This class method is principally intended to be called by the
        :meth:`YamlList.append_default` method, appending this instance to an
        existing parent list of similar instances.

        Parameters
        ----------
        yaml_list : YamlList
            Parent list instantiating this instance.
        '''

        pass


class YamlListItemTypedABC(YamlListItemABC):
    '''
    Abstract base class of all low-level YAML-backed typed list item subclasses,
    each instance of which encapsulates a YAML dictionary whose keys define the
    type and name of this item intended to be added to a :class:`YamlList`
    container.

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
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# CAUTION: To avoid diamond inheritance issues, neither this subclass nor any
# application-specific superclasses of this subclass (e.g., "YamlABC") should
# override any methods already defined by the "MutableSequence" superclass. This
# includes but is *NOT* limited to the public "MutableSequence" methods:
# append(), clear(), count(), extend(), index(), insert(), mro(), pop(),
# register(), remove(), and reverse().
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
class YamlList(YamlABC, MutableSequence):
    '''
    YAML-backed list subconfiguration both loaded from and saved back into a
    YAML-formatted file.

    Each item of this list is a :class:`YamlListItemABC` instance encapsulating
    the corresponding YAML-backed list item. Each such item may technically
    encapsulate any valid YAML type (e.g., :class:`int`, :class:`str`) but
    typically encapsulates a YAML dictionary of related key-value pairs (e.g.,
    animation settings, tissue profile).

    Attributes
    ----------
    _confs_wrap : list
        High-level list of each instance of the :attr:`_item_type` subclass
        encapsulating each dictionary in the low-level :attr:`_conf` list.
    _item_type : ClassType
        Subclass of the :class:`YamlListItemABC` abstract base class with
        which to instantiate each wrapper in the high-level :attr:`_confs_wrap`
        list encapsulating each dictionary in the low-level :attr:`_conf` list.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(self, item_type: ClassType) -> None:
        '''
        Initialize this low-level YAML-backed list.

        Attributes
        ----------
        item_type : ClassType
            Subclass of the :class:`YamlListItemABC` abstract base class
            with which to encapsulate each low-level YAML-backed list item in
            the passed ``confs`` list. For each such item, a new instance of
            this subclass is appended to the internal :attr:`_confs_wrap` list
            storing these instances.
        '''

        # Initialize our superclass.
        super().__init__()

        # If the passed type is *NOT* a high-level YAML-backed list item
        # subclass, raise an exception.
        classes.die_unless_subclass(
            subclass=item_type, superclass=YamlListItemABC)

        # Classify all passed parameters.
        self._item_type = item_type

        # Nullify all remaining parameters for safety.
        self._confs_wrap = None

    # ..................{ SUPERCLASS ~ YamlABC               }..................
    # Concrete methods provided by our YamlABC superclass.

    @type_check
    def load(self, conf: SequenceTypes) -> (
        'betse.lib.yaml.abc.yamllistabc.YamlList'):

        # Load our superclass with the passed sequence.
        super().load(conf=conf)

        # Initialize a high-level list of high-level wrappers.
        self._confs_wrap = []

        # For each low-level dictionary in this low-level sequence...
        for conf_yaml in conf:
            # High-level wrapper of the caller-specified type.
            conf_wrap = self._item_type()
            conf_wrap.load(conf=conf_yaml)

            # Wrap this dictionary with this wrapper.
            self._confs_wrap.append(conf_wrap)

        # Return this configuration for convenience.
        return self


    def unload(self) -> None:

        # Unload our superclass.
        super().unload()

        # Nullify this high-level list of high-level wrappers.
        self._confs_wrap = None

    # ..................{ SUPERCLASS ~ MutableSequence       }..................
    # Abstract methods required by our MutableSequence superclass.

    def __len__(self):
        return len(self._confs_wrap)


    @type_check
    def __delitem__(self, index: int) -> None:
        del self._conf[index]
        del self._confs_wrap[index]


    @type_check
    def __getitem__(self, index: int) -> YamlListItemABC:
        return self._confs_wrap[index]


    @type_check
    def __setitem__(self, index: int, value: YamlListItemABC) -> None:
        '''
        Set the low-level YAML-backed list item with the passed 0-based index to
        the passed such item.
        '''

        # If the passed object is *NOT* an instance of the caller-defined type,
        # raise an exception. This type is specified at initialization time and
        # hence cannot be type-checked via a method annotation.
        objects.die_unless_instance(obj=value, cls=self._item_type)

        # Set the high-level list item with this index to this object.
        self._confs_wrap[index] = value

        # Set the low-level list item with this index to this object's
        # underlying YAML-backed dictionary.
        self._conf[index] = value.conf


    @type_check
    def insert(self, index: int, value: YamlListItemABC) -> None:
        '''
        Insert the passed low-level YAML-backed list item immediately *before*
        the item with the passed 0-based index.
        '''

        # Raise an exception unless the passed object is an instance of the
        # desired API, which is specified at initialization time and hence
        # cannot be type checked above by a method annotation.
        objects.die_unless_instance(obj=value, cls=self._item_type)

        # Insert this object *BEFORE* the high-level list item with this index.
        self._confs_wrap.insert(index, value)

        # Insert this object's underlying YAML-backed dictionary *BEFORE* the
        # low-level list item with this index.
        self._conf.insert(index, value.conf)

    # ..................{ APPENDERS                          }..................
    def append_default(self) -> YamlListItemABC:
        '''
        Append a new low-level YAML-backed list item initialized to default
        values suitable for this list and return this item.

        Returns
        ----------
        YamlListItemABC
            Instance of the :attr:`_item_type` subclass appended to the
            high-level :attr:`_confs_wrap` list, encapsulating the new list item
            appended to the low-level :attr:`_conf` list.
        '''

        # New simulation configuration list item specific to this list.
        conf_wrap = self._item_type.make_default(self)

        # Append this list item to this list.
        self.append(conf_wrap)

        # Return this list item.
        return conf_wrap

    # ....................{ GETTERS                            }....................
    @type_check
    def get_item_name_unique(self, name_format: str) -> str:
        '''
        Identifier formatted according to the passed format specifier,
        guaranteed to be unique across all items of this list.

        By duck typing, each item of this list is required to define a ``name``
        instance variable whose value is the YAML-backed name of that item; if
        this is *not* the case, an exception is raised.

        This method then returns an identifier suitable for use by the caller as
        the value of the ``name`` instance variable for a newly created item of
        this list. Hence, callers tend to call this method from subclass
        :meth:`make_default` implementations.

        Parameters
        ----------
        name_format : str
            Format specifier containing a ``{}`` substring (e.g.,
            ``tissue ({})``), iteratively interpolated by this function with an
            arbitrary integer to produce this unique list item name.

        Returns
        ----------
        str
            Item name unique to this list matching this format.
        '''

        # Arbitrary unique identifier with which to uniquify (i.e., guarantee
        # the uniqueness of) the name of a new item in this list, defaulting to the
        # number of existing items in this list.
        yaml_list_item_id = len(self._confs_wrap)

        # Name of this tissue profile, unique in this list.
        yaml_list_item_name = None

        # While this name is *NOT* unique in this list, iteratively (re)search
        # this list until obtaining a unique name. This iteration is guaranteed
        # to (eventually) terminate successfully with a unique name.
        while True:
            yaml_list_item_name = name_format.format(yaml_list_item_id)

            # For each existing item of this list...
            for yaml_list_item_other in self._confs_wrap:
                # If this item already has this name, this name is non-unique.
                # In this case, increment this identifier, format a new name via
                # this identifier, and repeat this search.
                if yaml_list_item_name == yaml_list_item_other.name:
                    yaml_list_item_id += 1
                    break
            # If no item already has this name, this name is unique. In that
            # case, cease searching.
            else:
                break

        # Return this name.
        return yaml_list_item_name
