#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Abstract base classes of all YAML-backed simulation configuration subclasses as
well as functionality pertaining to such classes.
'''

# ....................{ IMPORTS                            }....................
from abc import ABCMeta, abstractmethod
from betse.lib.yaml import yamls
from betse.lib.yaml.yamlalias import yaml_alias
# FIXME: Ideally, submodules in the "betse.lib" subpackage should *NOT* import
# from submodules in the "betse.science" subpackage. To enforce this, the
# "SimPipeRunnerConfMixin" base class may need to be generalized into a generic
# mixin independent of simulations and then shifted into a new utility submodule
# (e.g., "betse.util.multi.piperun").
from betse.science.simulate.pipe.piperun import SimPipeRunnerConfMixin
from betse.util.io.log import logs
from betse.util.path import files, pathnames
from betse.util.type.cls import classes
from betse.util.type.obj import objects
from betse.util.type.types import (
    type_check,
    ClassType,
    MappingType,
    SequenceOrNoneTypes,
    StrOrNoneTypes,
)
from collections.abc import MutableSequence

# ....................{ SUPERCLASSES                       }....................
class YamlABC(object, metaclass=ABCMeta):
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
        Associate this high-level simulation configuration with the passed
        low-level dictionary.

        Parameters
        ----------
        conf : MappingType
            Low-level dictionary of related configuration settings both loaded
            from and savable back to the current YAML-formatted simulation
            configuration file.
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

# ....................{ SUPERCLASSES ~ file                }....................
class YamlFileABC(YamlABC):
    '''
    Abstract base class of all top-level simulation configuration subclasses,
    each directly backed by an low-level simulation configuration file in YAML
    format and hence both **serializable** (i.e., writable) to and
    **deserializable** (i.e., readable) from that file.

    Caveats
    ----------
    The ``conf_filename`` parameter accepted by most methods of this class
    (e.g., :meth:`read`, :meth:`write`) should be suffixed by a valid YAML
    filetype -- namely, either ``.yml`` or ``.yaml``. This class does *not*
    erforce this recommendation but does log non-fatal warnings where violated.

    External callers should ideally *never* access the public :meth:`conf`
    property returning a low-level dictionary containing this entire simulation
    configuration. The settings encapsulated by this dictionary are safely
    retrievable and modifiable by callers *only* via public :func:`yaml_alias`
    data descriptors leveraged by subclasses.

    Attributes
    ----------
    _conf_dirname : StrOrNoneTypes
        Absolute path of the directory containing the file whose filename is
        :attr:`conf_filename` and subdirectories with files referenced by this
        configuration file if such a file has been read (e.g., by a prior call
        to the :meth:`read` method) *or* ``None`` otherwise.
    _conf_filename : StrOrNoneTypes
        Absolute path of the low-level YAML-formatted simulation configuration
        file from which this object was most recently deserialized if such a
        file has been read (e.g., by a prior call to the :meth:`read` method)
        *or* ``None`` otherwise.
    '''

    # ..................{ MAKERS                             }..................
    @classmethod
    @type_check
    def make(cls, conf_filename: str) -> (
        'betse.lib.yaml.yamlabc.YamlFileABC'):
        '''
        Create return an instance of this subclass deserialized (i.e., read)
        from the passed YAML-formatted simulation configuration file.

        Parameters
        ----------
        conf_filename : str
            Absolute or relative path of the source file to be deserialized.
        '''

        self = cls()
        self.read(conf_filename)
        return self

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(self) -> None:
        '''
        Initialize this simulation configuration in the **closed state** (i.e.,
        associated with *no* low-level YAML-formatted simulation configuration
        file).
        '''

        # Initialize our superclass with the empty dictionary, which the
        # subsequent call to the read() method replaces with a non-empty
        # dictionary. While awkward, this approach avoids even *MORE* awkward
        # chicken-and-egg API issues.
        super().__init__(conf={})

        # Initialize this simulation configuration in the closed state. To avoid
        # extraneous logging, the body of the close() method is duplicated here.
        self._conf_dirname = None
        self._conf_filename = None

    # ..................{ PROPERTIES                         }..................
    @property
    def is_open(self) -> bool:
        '''
        ``True`` only if this simulation configuration is in the **open state**
        (i.e., associated with a low-level YAML-formatted simulation
        configuration file).

        If ``True``, *all* methods of this base class (e.g., :meth:`read`,
        :meth:`close`) are safely callable by callers; else, none are.
        '''

        return self._conf_filename is None

    # ..................{ CLOSERS                            }..................
    def close(self) -> None:
        '''
        Deassociate this high-level simulation configuration from its low-level
        YAML-formatted simulation configuration file if such a file has been
        read (e.g., by a prior call to the :meth:`read` method) *or* silently
        noop otherwise.
        '''

        # Log this operation.
        logs.log_info('Closing simulation configuration...')

        # Preserve the superclass contract that this variable be non-None.
        self._conf = {}

        # Nullify all instance variables for safety.
        self._conf_dirname = None
        self._conf_filename = None

    # ..................{ READERS                            }..................
    #FIXME: Validate the contents of this file (e.g., via "yamale").
    @type_check
    def read(self, conf_filename: str) -> None:
        '''
        Deserialize the passed YAML-formatted simulation configuration file into
        a low-level dictionary internally stored in this object, replacing the
        prior contents of this dictionary.

        Parameters
        ----------
        conf_filename : str
            Absolute or relative path of the source file to be deserialized.
        '''

        # Log this operation.
        logs.log_info(
            'Reading simulation configuration "%s"...',
            pathnames.get_basename(conf_filename))

        # Associate this object with this file.
        self._set_conf_filename(conf_filename)

        # Deserialize this file into this dictionary.
        self._conf = yamls.load(conf_filename)

    # ..................{ WRITERS                            }..................
    @type_check
    def overwrite(self) -> None:
        '''
        Serialize the low-level dictionary internally stored in this object to
        the current YAML-formatted simulation configuration file associated with
        this object, replacing the prior contents of this file.

        This method effectively implements the "Save" GUI metaphor.
        '''

        # Log this operation.
        logs.log_info('Overwriting simulation configuration...')

        # Delete this file (if found), preventing the subsequent write from
        # raising an otherwise ignorable exception.
        files.remove_if_found(self._conf_filename)

        # Resave this dictionary to this file.
        yamls.save(container=self._conf, filename=self._conf_filename)


    @type_check
    def write(self, conf_filename: str) -> None:
        '''
        Serialize the low-level dictionary internally stored in this object to
        the passed YAML-formatted simulation configuration file, replacing the
        prior contents of this file, *and* associate this object with this file.

        This method effectively implements the "Save As..." GUI metaphor.

        Parameters
        ----------
        conf_filename : str
            Absolute or relative path of the target file to be serialized.

        Raises
        ----------
        BetseFileException
            If this file already exists.
        '''

        # Log this operation.
        logs.log_info(
            'Writing simulation configuration "%s"...',
            pathnames.get_basename(conf_filename))

        # Validate this file *BEFORE* writing this file.
        files.die_if_file(conf_filename)

        # Associate this object with this file.
        self._set_conf_filename(conf_filename)

        # Save this dictionary to this file.
        yamls.save(container=self._conf, filename=self._conf_filename)

    # ..................{ PROPERTIES ~ read-only             }..................
    # Read-only properties, preventing callers from resetting these attributes.

    @property
    def conf_dirname(self) -> StrOrNoneTypes:
        '''
        Absolute path of the directory containing the file whose filename is
        :attr:`conf_filename` if such a file has been read (e.g., by a prior
        call to the :meth:`read` method) *or* ``None`` otherwise.
        '''

        return self._conf_dirname


    @property
    def conf_filename(self) -> StrOrNoneTypes:
        '''
        Absolute path of the low-level YAML-formatted simulation configuration
        file from which this object was most recently deserialized if such a
        file has been read (e.g., by a prior call to the :meth:`read` method)
        *or* ``None`` otherwise.
        '''

        return self._conf_filename

    # ..................{ SETTERS                            }..................
    @type_check
    def _set_conf_filename(self, conf_filename: str) -> None:
        '''
        Set the absolute path of the YAML-formatted file associated with this
        simulation configuration.

        Design
        ----------
        To prevent external callers from unsafely setting this path, this setter
        is intentionally implemented as an manual setter rather than a more
        preferable :meth:`conf_filename` property setter.
        '''

        # Unique absolute path of this file assigned *BEFORE* this file's
        # parent directory, ensuring the latter is non-empty.
        self._conf_filename = pathnames.canonicalize(conf_filename)

        # Unique absolute path of the parent directory of this file.
        self._conf_dirname = pathnames.get_dirname(self._conf_filename)

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
