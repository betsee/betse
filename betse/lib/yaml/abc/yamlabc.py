#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Abstract base classes of all YAML-backed simulation configuration subclasses as
well as functionality pertaining to such classes.
'''

# ....................{ IMPORTS                            }....................
from abc import ABCMeta
from betse.exceptions import BetseYamlException
from betse.lib.yaml import yamls
from betse.util.io.log import logs
from betse.util.path import dirs, pathnames
from betse.util.type.iterators import empty_iterator
from betse.util.type.types import (
    type_check,
    IterableTypes,
    MappingType,
    NoneType,
    StrOrNoneTypes,
)

# ....................{ SUPERCLASSES                       }....................
class YamlABC(object, metaclass=ABCMeta):
    '''
    Abstract base class of all YAML-backed configuration subclasses, each
    encapsulating a low-level YAML dictionary of related configuration settings
    (e.g., representing one tissue profile) both loaded from and savable back to
    a parent YAML-formatted configuration file.

    Attributes
    ----------
    _conf : MappingType
        Low-level dictionary of related configuration settings both loaded from
        and savable back to the parent YAML-formatted configuration file.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(self, conf: MappingType) -> None:
        '''
        Associate this high-level configuration with the passed low-level
        dictionary.

        Parameters
        ----------
        conf : MappingType
            Low-level dictionary of related configuration settings both loaded
            from and savable back to the parent YAML-formatted configuration
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
        savable back to the parent YAML-formatted configuration file.
        '''

        return self._conf


# Intended for use in callable type validation.
YamlABCOrNoneTypes = (YamlABC, NoneType)
'''
Tuple of both the YAML-backed configuration type *and* the type of the singleton
``None`` object.
'''

# ....................{ SUPERCLASSES ~ file                }....................
class YamlFileABC(YamlABC):
    '''
    Abstract base class of all top-level configuration subclasses, each directly
    backed by a low-level YAML-formatted configuration file in format and hence
    both **serializable** (i.e., writable) to and **deserializable** (i.e.,
    readable) from that file.

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
        Absolute path of the low-level YAML-formatted configuration file from
        which this object was most recently deserialized if any *or* ``None``
        otherwise.
    '''

    # ..................{ MAKERS                             }..................
    @classmethod
    @type_check
    def make(cls, conf_filename: str) -> (
        'betse.lib.yaml.abc.yamlabc.YamlFileABC'):
        '''
        Create return an instance of this subclass deserialized (i.e., read)
        from the passed YAML-formatted configuration file.

        Parameters
        ----------
        conf_filename : str
            Absolute or relative path of the source file to be deserialized.
        '''

        self = cls()
        self.load(conf_filename)
        return self

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(self) -> None:
        '''
        Initialize this file-backed configuration in the **unload state** (i.e.,
        associated with *no* low-level YAML-formatted configuration file).
        '''

        # Initialize our superclass with the empty dictionary, which the
        # subsequent call to the read() method replaces with a non-empty
        # dictionary. While awkward, this approach avoids even *MORE* awkward
        # chicken-and-egg API issues.
        super().__init__(conf={})

        # Initialize this simulation configuration in the unload state. To avoid
        # extraneous logging, the unload() method body is duplicated here.
        self._conf_dirname = None
        self._conf_filename = None

    # ..................{ PROPERTIES                         }..................
    @property
    def is_loaded(self) -> bool:
        '''
        ``True`` only if this simulation configuration is in the **loaded
        state** (i.e., deserialized from a low-level YAML-formatted simulation
        configuration file).

        If ``True``, *all* methods of this base class (e.g., :meth:`save`,
        :meth:`unload`,) are safely callable by callers; else, only the
        properties and the :meth:`read` method are.
        '''

        return self._conf_filename is not None

    # ..................{ LOADERS                            }..................
    #FIXME: Validate the contents of this file (e.g., via "yamale").
    @type_check
    def load(self, conf_filename: str) -> None:
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
            'Loading YAML file "%s"...', pathnames.get_basename(conf_filename))

        # Associate this object with this file.
        self._set_conf_filename(conf_filename)

        # Deserialize this file into this dictionary.
        self._conf = yamls.load(conf_filename)

    # ..................{ UNLOADERS                          }..................
    def unload(self) -> None:
        '''
        Deassociate this high-level simulation configuration from its low-level
        YAML-formatted simulation configuration file if such a file has been
        read (e.g., by a prior call to the :meth:`read` method) *or* silently
        noop otherwise.
        '''

        # Log this operation.
        logs.log_info('Closing YAML file...')

        # Preserve the superclass contract that this variable be non-None.
        self._conf = {}

        # Nullify all instance variables for safety.
        self._conf_dirname = None
        self._conf_filename = None

    # ..................{ SAVERS                             }..................
    @type_check
    def save(self, conf_filename: str) -> None:
        '''
        Serialize the low-level dictionary internally stored in this object to
        the passed YAML-formatted simulation configuration file, replacing the
        prior contents of this file, *and* associate this object with this file.

        This method effectively implements the "Save As..." GUI metaphor.

        Parameters
        ----------
        conf_filename : str
            Absolute or relative path of the target file to be serialized.
        '''

        # Log this operation.
        logs.log_info(
            'Saving YAML file "%s"...', pathnames.get_basename(conf_filename))

        # If no file to be saved has been read, raise an exception.
        self._die_unless_loaded()

        # Save this dictionary to this file.
        yamls.save(
            container=self._conf,
            filename=conf_filename,
            is_overwritable=True,
        )

        # Absolute paths of the parent directories containing the current file
        # and the passed file.
        src_dirname = self.conf_dirname
        trg_dirname = pathnames.get_dirname(conf_filename)

        # If these paths differ, all relative subdirectories internally
        # referenced and required by this file *MUST* be recursively copied from
        # the former to the latter.
        if src_dirname != trg_dirname:
            # For the absolute or relative path of each such subdirectory...
            for conf_subdirname in self._iter_conf_subdirnames():
                # If this path is absolute, silently ignore this subdirectory.
                if pathnames.is_absolute(conf_subdirname):
                    continue
                # Else, this path is relative and hence requires copying.

                # Absolute path of the old subdirectory.
                src_subdirname = pathnames.join(src_dirname, conf_subdirname)

                # Recursively copy from the old into the new subdirectory.
                dirs.copy_into_dir(
                    src_dirname=src_subdirname,
                    trg_dirname=trg_dirname,
                    is_overwritable=True,
                )

        # Associate this object with this file *AFTER* successfully copying to
        # this file and all external paths required by this file.
        self._set_conf_filename(conf_filename)


    @type_check
    def save_inplace(self) -> None:
        '''
        Serialize the low-level dictionary internally stored in this object to
        the current YAML-formatted simulation configuration file associated with
        this object, replacing the prior contents of this file.

        This method effectively implements the "Save" GUI metaphor.
        '''

        # Log this operation.
        logs.log_info('Overwriting YAML file...')

        # If no file to be saved has been read, raise an exception.
        self._die_unless_loaded()

        # Resave this dictionary to this file.
        yamls.save(
            container=self._conf,
            filename=self._conf_filename,
            is_overwritable=True,
        )

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

    # ..................{ EXCEPTIONS                         }..................
    def _die_unless_loaded(self) -> None:
        '''
        Raise an exception unless this file-backed configuration is currently in
        the **read state** (i.e., associated with a low-level YAML-formatted
        configuration file).
        '''

        if not self.is_loaded:
            raise BetseYamlException('No YAML file open.')

    # ..................{ SETTERS                            }..................
    @type_check
    def _set_conf_filename(self, conf_filename: str) -> None:
        '''
        Set the absolute path of the YAML-formatted file associated with this
        configuration.

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

    # ..................{ SUBCLASS ~ optional                }..................
    # Methods intended to be optionally overriden by subclasses.

    def _iter_conf_subdirnames(self) -> IterableTypes:
        '''
        Generator yielding the pathname of each requisite direct subdirectory of
        the parent directory of the YAML-formatted file currently associated
        with this configuration, where "requisite" means internally referenced
        and hence required by this file.

        For safety, this generator excludes all direct subdirectories of this
        parent directory that are *not* internally referenced by this file.

        By default, this superclass method returns the empty generator.

        Yields
        ----------
        str
            Absolute or relative pathname of each such subdirectory.
        '''

        # If no file has been read, raise an exception.
        self._die_unless_loaded()

        # Default to the empty iterator.
        return empty_iterator()
