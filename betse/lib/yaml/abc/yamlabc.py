#!/usr/bin/env python3
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
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
    MappingOrSequenceTypes,
    NoneType,
)

# ....................{ SUPERCLASSES                       }....................
class YamlABC(object, metaclass=ABCMeta):
    '''
    Abstract base class of all YAML-backed in-memory configuration subclasses,
    each encapsulating a low-level container of related configuration settings
    (e.g., representing one tissue profile) both loaded from and saved back to a
    parent YAML-formatted configuration file.

    Attributes
    ----------
    _conf : MappingOrSequenceOrNoneTypes
        Low-level container of related configuration settings loaded from and
        saved back to a YAML-formatted configuration file if the :meth:`load`
        method has been called *or* ``None`` otherwise.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(self) -> None:
        '''
        Initialize this YAML-backed configuration.
        '''

        self._conf = None

    # ..................{ PROPERTIES ~ read-only             }..................
    # Read-only properties, preventing callers from resetting these attributes.

    @property
    def is_loaded(self) -> bool:
        '''
        ``True`` only if this configuration is in the **loaded state** (i.e.,
        the :meth:`load` method has been called more recently than the
        :meth:`unload` method has been called).

        If ``True``, the :meth:`conf` property is safely accessible by callers.
        '''

        return self._conf is not None


    @property
    def conf(self) -> MappingOrSequenceTypes:
        '''
        Low-level container of related configuration settings loaded from and
        saved back to a YAML-formatted configuration file if the :meth:`load`
        method has been called more recently than the :meth:`unload` method has
        been called *or* raise an exception otherwise.

        Raises
        ----------
        BetseYamlException
             If the :meth:`load` method has yet be called.
        '''

        # If this property has yet to be set, raise an exception.
        self.die_unless_loaded()

        # Else, this property is a container. Return this container.
        return self._conf

    # ..................{ EXCEPTIONS                         }..................
    def die_unless_loaded(self) -> None:
        '''
        Raise an exception unless this configuration is in the loaded state.

        See Also
        ----------
        :meth:`is_loaded`
            Further details.
        '''

        if not self.is_loaded:
            raise BetseYamlException(
                'YAML configuration not loaded '
                '(i.e., load() method not called yet).')

    # ..................{ LOADERS                            }..................
    @type_check
    def load(self, conf: MappingOrSequenceTypes) -> (
        'betse.lib.yaml.abc.yamlabc.YamlABC'):
        '''
        Associate this configuration with the passed YAML-backed container.

        Parameters
        ----------
        conf : MappingOrSequenceTypes
            Low-level container of related configuration settings loaded from
            and saved back to a YAML-formatted configuration file.

        Returns
        ----------
        YamlABC
            This configuration, returned for convenience chaining in callers.
        '''

        # Classify this container.
        self._conf = conf

        # Return this configuration for convenience.
        return self


    def unload(self) -> None:
        '''
        Deassociate this configuration from its YAML-backed container if any
        *or* reduce to a noop otherwise.

        This method serves as a low-level safety mechanism ensuring that the
        next access of the :meth:`conf` property will raise a human-readable
        exception.
        '''

        self._conf = None


# Intended for use in callable type validation.
YamlABCOrNoneTypes = (YamlABC, NoneType)
'''
Tuple of both the YAML-backed configuration type *and* the type of the singleton
``None`` object.
'''

# ....................{ SUPERCLASSES ~ file                }....................
class YamlFileABC(YamlABC):
    '''
    Abstract base class of all YAML-backed in-memory and on-disk configuration
    subclasses, encapsulating a low-level container of *all* configuration
    settings both loaded from and saved back to a YAML-formatted configuration
    file.

    Caveats
    ----------
    The ``conf_filename`` parameter accepted by most methods of this class
    (e.g., :meth:`load`, :meth:`write`) should be suffixed by a valid YAML
    filetype -- namely, either ``.yml`` or ``.yaml``. This class does *not*
    erforce this recommendation but does log non-fatal warnings where violated.

    External callers should *never* access the low-level :meth:`conf` property,
    as doing so violates type safety and value constraints. The configuration
    settings persisted by this property are safely retrievable and modifiable
    *only* via :func:`yaml_alias` data descriptors defined by subclasses.

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

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(self) -> None:
        '''
        Initialize this file-backed configuration in the **unload state** (i.e.,
        associated with *no* low-level YAML-formatted configuration file).
        '''

        # Initialize our superclass.
        super().__init__()

        # Initialize this simulation configuration in the unload state. To avoid
        # extraneous logging, the unload() method body is duplicated here.
        self._conf_dirname = None
        self._conf_filename = None

    # ..................{ PROPERTIES ~ read-only             }..................
    # Read-only properties, preventing callers from resetting these attributes.

    @property
    def conf_filename(self) -> str:
        '''
        Absolute path of the low-level YAML-formatted file from which this file
        was deserialized if such a file has been loaded by a prior call to the
        :meth:`load` method *or* raise an exception otherwise.

        Raises
        ----------
        BetseYamlException
             If the :meth:`load` method has yet to be called.
        '''

        # If this property has yet to be set, raise an exception.
        self.die_unless_loaded()

        # Return this property.
        return self._conf_filename


    @property
    def conf_dirname(self) -> str:
        '''
        Absolute path of the directory containing the file with filename
        :attr:`conf_filename` if such a file has been loaded by a prior call to
        the :meth:`load` method *or* raise an exception otherwise.

        Raises
        ----------
        BetseYamlException
             If the :meth:`load` method has yet to be called.
        '''

        # If this property has yet to be set, raise an exception.
        self.die_unless_loaded()

        # Return this property.
        return self._conf_dirname

    # ..................{ LOADERS                            }..................
    @type_check
    def load(self, conf_filename: str) -> (
        'betse.lib.yaml.abc.yamlabc.YamlFileABC'):
        '''
        Deserialize the passed YAML-formatted simulation configuration file into
        a low-level dictionary internally stored in this object, replacing the
        prior contents of this dictionary.

        Parameters
        ----------
        conf_filename : str
            Absolute or relative path of the source file to be deserialized.

        Returns
        ----------
        YamlFileABC
            This configuration, returned for convenience chaining in callers.
        '''

        # Log this operation.
        logs.log_debug(
            'Loading YAML file "%s"...', pathnames.get_basename(conf_filename))

        # Low-level dictionary deserialized from this file.
        conf = yamls.load(filename=conf_filename)

        # Load this dictionary into our superclass.
        super().load(conf=conf)

        # Associate this object with this file *AFTER* successfully
        # deserializing this file.
        self._set_conf_filename(conf_filename)

        # Return this configuration for convenience.
        return self


    def unload(self) -> None:

        # If this file is actually loaded, log this operation.
        if self.is_loaded:
            logs.log_debug(
                'Closing YAML file "%s"...',
                pathnames.get_basename(self.conf_filename))

        # Unload our superclass.
        super().unload()

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
        logs.log_debug(
            'Saving YAML file "%s"...', pathnames.get_basename(conf_filename))

        # Save this dictionary to this file.
        yamls.save(
            container=self.conf,
            filename=conf_filename,
            is_overwritable=True,
        )

        # Absolute paths of the parent directories containing the current file
        # and the passed file, canonicalized to permit comparison below.
        src_dirname = pathnames.canonicalize(self.conf_dirname)
        trg_dirname = pathnames.canonicalize(
            pathnames.get_dirname(conf_filename))

        # If these paths differ, recursively copy all relative subdirectories
        # internally referenced and hence required by this file.
        if src_dirname != trg_dirname:
            # For relative pathname of each such subdirectory...
            for conf_subdirname in self._iter_conf_subdirnames():
                # Absolute pathname of the old subdirectory.
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
        logs.log_debug(
            'Overwriting YAML file "%s"...',
            pathnames.get_basename(self.conf_filename))

        # Resave this dictionary to this file.
        yamls.save(
            container=self.conf,
            filename=self.conf_filename,
            is_overwritable=True,
        )

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
        self.die_unless_loaded()

        # Default to the empty iterator.
        return empty_iterator()
