#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Abstract base classes of all YAML-backed configuration subclasses encapsulating
on-disk lists and list items.
'''

# ....................{ IMPORTS                           }....................
from abc import abstractmethod
from betse.lib.yaml.abc.yamlabc import YamlABC
from betse.lib.yaml.abc.yamlmixin import YamlTypedBooledMixin
from betse.util.io.log import logs
from betse.util.type.cls import classes
from betse.util.type.obj import objects
from betse.util.type.text.string import strs
from betse.util.type.types import type_check, ClassType, SequenceTypes
from collections.abc import MutableSequence

# ....................{ SUPERCLASSES ~ list item          }....................
class YamlListItemABC(YamlABC):
    '''
    Abstract base class of all **YAML-backed list item** (i.e., configuration
    intended to be added to a parent :class:`YamlList` container) subclasses.

    Each instance of this superclass may technically encapsulate any valid YAML
    type (e.g., :class:`int`, :class:`str`) but typically encapsulates a YAML
    dictionary of related key-value pairs (e.g., settings configuring a single
    animation or tissue profile).
    '''

    # ..................{ MAKERS                            }..................
    @classmethod
    def make_list(cls, *args, **kwargs) -> (
        'betse.lib.yaml.abc.yamllistabc.YamlList'):
        '''
        Create and return a YAML-backed list containing only items of this
        subclass type, both loaded from and savable back to a YAML-formatted
        file.

        Parameters
        ----------
        All parameters are passed as is to the :meth:`YamlList.__init__`
        method.
        '''

        return YamlList(*args, item_type=cls, **kwargs)

    # ..................{ MAKERS ~ abstract                 }..................
    # Subclasses are required to implement the following abstract class
    # methods.

    @classmethod
    @abstractmethod
    def make_default(
        cls, yaml_list: 'betse.lib.yaml.abc.yamllistabc.YamlList',
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


class YamlListItemTypedBooledABC(YamlTypedBooledMixin, YamlListItemABC):
    '''
    Abstract base class of all **YAML-backed typed booled list item** (i.e.,
    YAML-backed list item whose value is a YAML dictionary defining a top-level
    key ``type`` whose value is a machine-readable string identifying that
    item's type *and* a top-level key ``enabled`` whose value is a boolean
    specifying whether that item is enabled or disabled) subclasses.
    '''

    pass

# ....................{ SUPERCLASSES ~ list               }....................
class YamlList(YamlABC, MutableSequence):
    '''
    YAML-backed list subconfiguration both loaded from and saved back into a
    YAML-formatted file.

    Each item of this list is a :class:`YamlListItemABC` instance encapsulating
    the corresponding YAML-backed list item. Each such item may technically
    encapsulate any valid YAML type (e.g., :class:`int`, :class:`str`) but
    typically encapsulates a YAML dictionary of related key-value pairs (e.g.,
    animation settings, tissue profile).

    Design
    ----------
    This list conforms to the :class:`MutableSequence` API and hence defines
    all public methods mandated by that API, including the familiar
    :meth:`append`, :meth:`clear`, :meth:`count`, :meth:`extend`,
    :meth:`index`, :meth:`insert`, :meth:`pop`, :meth:`register`,
    :meth:`remove`, and :meth:`reverse` methods.

    Attributes
    ----------
    _conf : SequenceTypes
        Low-level sequence of related configuration settings loaded from and
        saved back to a YAML-formatted configuration file if the :meth:`load`
        method has been called *or* ``None`` otherwise.
    _confs_wrap : list
        High-level list of each instance of the :attr:`_item_type` subclass
        encapsulating each dictionary in the low-level :attr:`_conf` sequence.
    _item_type : ClassType
        Subclass of the :class:`YamlListItemABC` abstract base class with
        which to instantiate each wrapper in the high-level :attr:`_confs_wrap`
        list encapsulating each dictionary in the low-level :attr:`_conf` list.
    '''

    # ..................{ INITIALIZERS                      }..................
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

    # ..................{ SUPERCLASS ~ YamlABC              }..................
    # Concrete methods provided by our YamlABC superclass.

    @type_check
    def load(self, conf: SequenceTypes) -> None:

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


    def unload(self) -> None:

        # Unload our superclass.
        super().unload()

        # Nullify this high-level list of high-level wrappers.
        self._confs_wrap = None

    # ..................{ SUPERCLASS ~ mutable : abstract   }..................
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
        Set the low-level YAML-backed list item with the passed 0-based index
        to the passed such item.
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

    # ..................{ SUPERCLASS ~ mutable : concrete   }..................
    # Concrete methods defined by our MutableSequence superclass and overridden
    # here for dubious reasons (e.g., premature optimization).

    def clear(self) -> None:
        '''
        Remove all YAML-backed list items from this list, reducing this list to
        to the empty list.
        '''

        # Clear both the low- and high-level sequences wrapped by this object.
        self._conf.clear()
        self._confs_wrap.clear()

    # ..................{ APPENDERS                         }..................
    def append_default(self) -> YamlListItemABC:
        '''
        Append a new low-level YAML-backed list item initialized to default
        values suitable for this list *and* return this item.

        Returns
        ----------
        YamlListItemABC
            Instance of the :attr:`_item_type` subclass appended to the
            high-level :attr:`_confs_wrap` list, encapsulating the new list
            item appended to the low-level :attr:`_conf` list.
        '''

        # New simulation configuration list item specific to this list.
        conf_wrap = self._item_type.make_default(yaml_list=self)

        # Append this list item to this list.
        self.append(conf_wrap)

        # Return this list item.
        return conf_wrap

    # ..................{ GETTERS                           }..................
    #FIXME: Generalize this method into a new
    #betse.util.type.iterable.iterables.uniquify_item_attr_value() function
    #accepting a new "iterable: IterableTypes" first parameter for generality.
    #For convenience, retain this method rewritten in terms of that newly
    #defined function.
    @type_check
    def uniquify_item_attr_value(
        self, attr_name: str, attr_value_format: str) -> str:
        '''
        Create and return a machine-readable string guaranteed to be unique
        across all items of this list when set as the value of the settable
        attribute (e.g., YAML alias) with the passed name of any such item,
        including an item that has yet to be created.

        This method requires each item of this list to declare an attribute
        with the passed name whose value is an arbitrary string. This method
        then returns a string suitable for use by the caller as the value of
        that attribute for a newly created item of this list. Hence, callers
        typically call this method from subclass :meth:`make_default`
        implementations.

        Parameters
        ----------
        attr_name : str
            Name of the attribute declared by all items of this list whose
            string value is to be uniquified.
        attr_value_format : str
            Format specifier containing a ``{}`` substring (e.g.,
            ``tissue ({})``), iteratively interpolated by this function with an
            arbitrary integer to produce this unique list item attribute value.

        Returns
        ----------
        str
            Value of this attribute guaranteed to both match this format *and*
            be unique across all items of this list.

        Raises
        ----------
        BetseStrException
            If this format specifier contains no ``{}`` substring.
        '''

        # Avoid circular import dependencies.
        from betse.util.type.obj import objects

        # Log this formatting.
        logs.log_debug(
            'Uniquifying YAML list item key "%s" with template "%s"...',
            attr_name, attr_value_format)

        #FIXME: Ideally, we would also raise exceptions if this format
        #specifier contains two or more ``{}`` substrings. Sadly, there appears
        #no trivial and/or efficient means of doing so. While we suppose we
        #*COULD* define a new
        #betse.util.type.text.string.strs.get_substrs_count() function doing
        #so, the lack of any benefit to doing so hardly seems worth it. Note
        #also that genuinely testing this is effectively infeasible due to
        #false positives (e.g., the string "{}{}" should raise exceptions but
        #the string "{}{{}}" should *NOT*).

        # If this specifier contains no "{}" substring, raise an exception.
        #
        # Note that this simplistic logic fails to account for "{{" and "}}"
        # escaping and hence *COULD* fail to raise exceptions when passed
        # worst-case format specifiers, but that we mostly do not care.
        strs.die_unless_substr(text=attr_value_format, substr='{}')

        # Arbitrary unique identifier with which to uniquify (i.e., guarantee
        # the uniqueness of) an arbitrary string attribute of a new item in
        # this list, defaulting to the number of existing items in this list.
        yaml_list_item_id = len(self._confs_wrap)

        # Arbitrary string attribute of a new item in this list.
        yaml_list_item_attr = None

        # While this attribute is *NOT* unique across this list, iteratively
        # (re)search this list until obtaining a unique attribute. This
        # iteration is guaranteed to terminate successfully.
        while True:
            yaml_list_item_attr = attr_value_format.format(yaml_list_item_id)

            # For each existing item of this list...
            for yaml_list_item_other in self._confs_wrap:
                # String value of the attribute with this name defined by this
                # item if such an attribute exists *OR* raise an exception.
                yaml_list_item_other_attr = objects.get_attr(
                    obj=yaml_list_item_other,
                    attr_name=attr_name,
                    attr_type=str,
                )

                # If this item already has this attribute, this attribute is
                # non-unique. In this case, increment this identifier, format a
                # new attribute via this identifier, and continue this search.
                if yaml_list_item_attr == yaml_list_item_other_attr:
                    yaml_list_item_id += 1
                    break
            # If no item already has this attribute, this attribute is unique.
            # In this case, cease searching.
            else:
                break

        # Return this attribute.
        return yaml_list_item_attr
