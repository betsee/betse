#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Abstract base classes of all YAML-backed file wrapper subclasses.
'''

# ....................{ IMPORTS                           }....................
from betse.lib.yaml import yamls
from betse.lib.yaml.abc.yamlabc import YamlABC
from betse.util.io.log import logs
from betse.util.path import dirs, pathnames
from betse.util.path.dirs import DirOverwritePolicy
from betse.util.type.descriptor.descs import abstractclassproperty_readonly
from betse.util.type.iterable.iterators import empty_iterator
from betse.util.type.types import type_check, IterableTypes

# ....................{ SUPERCLASSES                      }....................
class YamlFileABC(YamlABC):
    '''
    Abstract base class of all **YAML-backed file wrapper** (i.e., high-level
    object wrapping a low-level mapping or sequence of YAML-backed
    configuration settings both loadable from and savable back to a
    YAML-formatted file) subclasses.

    Caveats
    ----------
    The ``conf_filename`` parameter accepted by most methods of this class
    (e.g., :meth:`load`, :meth:`save`) should be suffixed by a valid YAML
    filetype -- namely, either ``.yml`` or ``.yaml``. This class does *not*
    erforce this recommendation but does log non-fatal warnings where violated.

    External callers should *never* access the low-level :meth:`conf` property,
    as doing so violates type safety and value constraints. The configuration
    settings persisted by this property are safely retrievable and modifiable
    *only* via :func:`yaml_alias` data descriptors defined by subclasses.

    Attributes
    ----------
    _conf_basename : StrOrNoneTypes
        Basename of :attr:`_conf_filename` if non-``None`` *or* ``None``
        otherwise.
    _conf_dirname : StrOrNoneTypes
        Absolute dirname of :attr:`_conf_filename` if non-``None`` *or*
        ``None`` otherwise.
    _conf_filename : StrOrNoneTypes
        Absolute filename of the low-level YAML-formatted configuration file
        from which this object was most recently deserialized if any *or*
        ``None`` otherwise.
    '''

    # ..................{ MAKERS                            }..................
    @classmethod
    @type_check
    def make(cls, conf_filename: str, *args, **kwargs) -> (
        'betse.lib.yaml.abc.yamlfileabc.YamlFileABC'):
        '''
        Create and return a YAML file wrapper of this subclass type,
        deserialized from the passed YAML-formatted file into a low-level
        mapping or sequence internally persisted in this wrapper.

        Parameters
        ----------
        conf_filename : str
            Absolute or relative filename of the source file to be
            deserialized.

        All other parameters are passed as is to the :meth:`__init__` method.

        Returns
        ----------
        YamlFileABC
            This configuration, returned for convenience chaining in callers.
        '''

        # Instance of this subclass, passed all passed variadic parameters.
        yaml_file_conf = cls(*args, **kwargs)

        # Deserialize the passed file into this instance.
        yaml_file_conf.load(conf_filename)

        # Return this instance.
        return yaml_file_conf

    # ..................{ INITIALIZERS                      }..................
    def __init__(self) -> None:
        '''
        Initialize this YAML file wrapper in the **unload state** (i.e.,
        associated with *no* low-level YAML-formatted configuration file).
        '''

        # Initialize our superclass.
        super().__init__()

        # Initialize this simulation configuration in the unload state. To
        # avoid extraneous logging, the unload() method body is repeated here.
        self._conf_basename = None
        self._conf_dirname = None
        self._conf_filename = None

    # ..................{ PROPERTIES ~ read-only            }..................
    # Read-only properties, preventing callers from resetting these attributes.

    @property
    def conf_basename(self) -> str:
        '''
        Basename of :attr:`conf_filename` if a prior call to the :meth:`load`
        method successfully loaded such a file *or* raise an exception
        otherwise.

        Raises
        ----------
        BetseYamlException
             If the :meth:`load` method has yet to be called.
        '''

        # If this property has yet to be set, raise an exception.
        self.die_unless_loaded()

        # Return this property.
        return self._conf_basename


    @property
    def conf_filename(self) -> str:
        '''
        Absolute filename of the low-level YAML-formatted file from which this
        file was deserialized if a prior call to the :meth:`load` method
        successfully loaded such a file *or* raise an exception otherwise.

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
        Absolute direname of :attr:`conf_filename` if a prior call to the
        :meth:`load` method successfully loaded such a file *or* raise an
        exception otherwise.

        Raises
        ----------
        BetseYamlException
             If the :meth:`load` method has yet to be called.
        '''

        # If this property has yet to be set, raise an exception.
        self.die_unless_loaded()

        # Return this property.
        return self._conf_dirname

    # ..................{ LOADERS                           }..................
    @type_check
    def load(self, conf_filename: str, **kwargs) -> None:
        '''
        Deserialize (i.e., load, read) the passed YAML-formatted file into a
        low-level mapping or sequence internally persisted in this wrapper.

        This method may be safely called multiple times with different files.
        This method implicitly calls the :meth:`unload` method if this is
        *not* the first call to this method, thus automatically deassociating
        this wrapper from all previously loaded files as needed.

        Parameters
        ----------
        conf_filename : str
            Absolute or relative filename of the source file to be
            deserialized.

        All remaining keyword arguments are passed as is to the
        :func:`betse.lib.yaml.yamls.load` function.
        '''

        # If a file is already loaded, unload this file for safety. While doing
        # so effectively reduces to a noop with respect to this superclass,
        # subclass implementations of the unload() method may elect to perform
        # critical functionality that *MUST* be respected here.
        if self.is_loaded:
            self.unload()

        # Log this load.
        logs.log_debug('Loading YAML file: %s', conf_filename)

        # Low-level dictionary deserialized from this file.
        conf = yamls.load(filename=conf_filename, **kwargs)

        # Load this dictionary into our superclass.
        super().load(conf=conf)

        # Associate this object with this file *AFTER* successfully
        # deserializing this file.
        self._set_conf_filename(conf_filename)


    def unload(self) -> None:

        # If this file is actually loaded, log this operation.
        if self.is_loaded:
            logs.log_debug('Closing YAML file: %s', self._conf_filename)

        # Unload our superclass.
        super().unload()

        # Nullify all instance variables for safety.
        self._conf_basename = None
        self._conf_dirname = None
        self._conf_filename = None

    # ..................{ SAVERS                            }..................
    #FIXME: For brevity, rename the following parameters:
    #
    #* "is_conf_file_overwritable" to "is_file_overwritable".
    #* "conf_subdir_overwrite_policy" to "subdir_overwrite_policy".
    @type_check
    def save(
        self,

        # Mandatory parameters.
        conf_filename: str,

        # Optional parameters.
        is_conf_file_overwritable: bool = False,
        conf_subdir_overwrite_policy: DirOverwritePolicy = (
            DirOverwritePolicy.SKIP_WITH_WARNING),
    ) -> None:
        '''
        Serialize (i.e., save, write) the low-level mapping or sequence
        internally persisted in this wrapper to the YAML-formatted file with
        the passed filename, copying *all* external resources internally
        referenced by this mapping or sequence into this file's directory.

        This method effectively implements the "Save As..." GUI metaphor.
        Specifically, this method (in order):

        #. Serializes this mapping or sequence to the file with this filename,
           optionally overwriting the existing contents of this file depending
           on the passed ``is_conf_file_overwritable`` parameter.
        #. Recursively copies all relative subdirectories internally referenced
           (and hence required) by this file from the directory of the current
           file associated with this wrapper into the directory of the passed
           file, optionally overwriting the existing contents of these
           subdirectories depending on the passed
           ``conf_subdir_overwrite_policy`` parameter.
        #. Associates this wrapper with this filename.

        Parameters
        ----------
        conf_filename : str
            Absolute or relative filename of the target file to be serialized.
        is_conf_file_overwritable : optional[bool]
            If this target file already exists *and* this boolean is:

            * ``True``, this target file is silently overwritten.
            * ``False``, an exception is raised.

            Defaults to ``False``.
        conf_subdir_overwrite_policy : DirOverwritePolicy
            **Subdirectory overwrite policy** (i.e., strategy for handling
            existing subdirectories to be overwritten by this save), where the
            subdirectories in question are all subdirectories of the directory
            of this target file. Defaults to
            :attr:`DirOverwritePolicy.SKIP_WITH_WARNING`, ignoring each target
            subdirectory that already exists with a non-fatal warning.

        Raises
        ----------
        BetseDirException
            If the passed ``conf_subdir_overwrite_policy`` parameter is
            :attr:`DirOverwritePolicy.HALT_WITH_EXCEPTION` *and* one or more
            subdirectories of the target directory already exist that are also
            subdirectories of the source directory.
        BetseFileException
            If the passed ``is_conf_file_overwritable`` parameter is ``False``
            *and* this target file already exists.
        '''

        # Log this save.
        logs.log_debug(
            'Saving YAML file: %s -> %s',
            self._conf_filename, conf_filename)

        # Serialize this mapping or sequence to this file.
        #
        # Note that the die_unless_loaded() method need *NOT* be called here,
        # as the "conf" property implicitly does so on our behalf.
        yamls.save(
            container=self.conf,
            filename=conf_filename,
            is_overwritable=is_conf_file_overwritable,
        )

        # Absolute dirnames of the directories containing the current file and
        # the passed file, canonicalized to permit comparison below.
        src_dirname = pathnames.canonicalize(self.conf_dirname)
        trg_dirname = pathnames.canonicalize(
            pathnames.get_dirname(conf_filename))

        # If these directories differ, recursively copy all relative
        # subdirectories internally referenced and hence required by this file.
        if src_dirname != trg_dirname:
            # For the basename of each such subdirectory...
            #
            # Note that the ideal solution of recursively copying this source
            # directory into the directory of this target file (e.g., via
            # "dirs.copy(src_dirname, pathnames.get_dirname(conf_filename))")
            # fails for the following subtle reasons:
            #
            # * This target directory may be already exist, which dirs.copy()
            #   prohibits even when the directory is empty.
            # * This target configuration file basename may differ from that of
            #   this source configuration file, necessitating a subsequent call
            #   to file.move().
            for conf_subdirname in self._iter_conf_subdir_basenames():
                # Absolute dirname of the source subdirectory.
                src_subdirname = pathnames.join(src_dirname, conf_subdirname)

                # Recursively copy from the old into the new subdirectory.
                dirs.copy_dir_into_dir(
                    src_dirname=src_subdirname,
                    trg_dirname=trg_dirname,
                    overwrite_policy=conf_subdir_overwrite_policy,

                    # Ignore all empty ".gitignore" files in all subdirectories
                    # of this source subdirectory. These files serve only as
                    # placeholders instructing Git to track otherwise empty
                    # subdirectories. Preserving such files only invites end
                    # user confusion.
                    ignore_basename_globs=('.gitignore',),
                )

        # Associate this object with this file *AFTER* successfully copying to
        # this file and all external paths required by this file.
        self._set_conf_filename(conf_filename)


    def save_inplace(self) -> None:
        '''
        Serialize the low-level dictionary internally stored in this object to
        the current YAML-formatted simulation configuration file associated
        with this object, replacing the prior contents of this file.

        This method effectively implements the "Save" GUI metaphor.
        '''

        # Log this save.
        logs.log_debug('Resaving YAML file: %s', self._conf_filename)

        # Resave this dictionary to this file.
        yamls.save(
            container=self.conf,
            filename=self.conf_filename,
            is_overwritable=True,
        )

    # ..................{ COPIERS                           }..................
    def copy(
        self,
        src_conf_filename: str,
        trg_conf_filename: str,
        **kwargs
    ) -> None:
        '''
        Copy the source to target YAML-formatted file with the passed filenames
        and *all* external resources from the directory of the former into that
        of the latter.

        Specifically, this method (in order):

        #. Deassociates this wrapper from any previously loaded file.
        #. Loads the passed source file into this wrapper.
        #. Saves this source file to the passed target file.
        #. Associates this wrapper with the passed target file.

        Note that this method is currently syntactic sugar for the following
        equivalent operations, where ``p`` is the current object:

        #. Calling ``p.load(src_conf_filename)``.
        #. Calling ``p.save(trg_conf_filename, **kwargs)``.

        Parameters
        ----------
        src_conf_filename : str
            Absolute or relative filename of the source YAML-formatted file to
            safely copy from.
        trg_conf_filename : str
            Absolute or relative filename of the target YAML-formatted file to
            safely copy from.

        All remaining keyword arguments are passed to the :meth:`save` method.
        '''

        # Load this source file into this wrapper, implicitly deassociating
        # this wrapper from any previously loaded file.
        self.load(src_conf_filename)

        # Save this source file to this target file.
        self.save(trg_conf_filename, **kwargs)

    # ..................{ SETTERS                           }..................
    @type_check
    def _set_conf_filename(self, conf_filename: str) -> None:
        '''
        Set the absolute path of the YAML-formatted file associated with this
        configuration.

        Design
        ----------
        To prevent external callers from unsafely setting this path, this
        setter is intentionally implemented as an manual setter rather than a
        more preferable :meth:`conf_filename` property setter.
        '''

        # Unique absolute filename of this file assigned *BEFORE* this file's
        # parent directory, ensuring the latter is non-empty.
        self._conf_filename = pathnames.canonicalize(conf_filename)

        # Basename of this file.
        self._conf_basename = pathnames.get_basename(self._conf_filename)

        # Unique absolute dirname of the parent directory of this file.
        self._conf_dirname = pathnames.get_dirname(self._conf_filename)

    # ..................{ SUBCLASS ~ optional               }..................
    # Methods intended to be optionally overriden by subclasses.

    def _iter_conf_subdir_basenames(self) -> IterableTypes:
        '''
        Generator iteratively yielding the basename of each **requisite
        subdirectory** (i.e., subdirectory internally referenced and hence
        required by the YAML-formatted file currently associated with this
        configuration) of the directory containing this file.

        For safety, this generator excludes all direct subdirectories of this
        parent directory that are *not* internally referenced by this file.

        By default, this superclass method returns the empty generator.

        Yields
        ----------
        str
            Basename of each such subdirectory.
        '''

        # If no file has been read, raise an exception.
        self.die_unless_loaded()

        # Default to the empty iterator.
        return empty_iterator()

# ....................{ SUPERCLASSES ~ default            }....................
#FIXME: Is this genuinely required anymore? It would probably be saner to
#simply refactor all remaining calls to the copy_default() method to call the
#copy() method instead with the default simulation configuration file.
#
#The current approach is a bit too "magical." In particular, the "Parameter"
#subclass' implementation of the "conf_default_filename" property assumes
#the application metadata singleton to have already been initialized, inviting
#complications whenever either that property or the copy_default() method are
#called. Since explicit is better than implicit, our usage of properties here
#is inadvisable at best and bug-prone at worst.
class YamlFileDefaultABC(YamlFileABC):
    '''
    Abstract base class of all **YAML-backed defaultable file wrapper** (i.e.,
    high-level object wrapping a low-level mapping or sequence of YAML-backed
    configuration settings both loadable from and savable back to a
    YAML-formatted file, optionally defaulting to a subclass-specific default
    YAML-formatted file) subclasses.

    This superclass augments the :class:`YamlFileABC` superclass with support
    for loading from a subclass-specific default YAML-formatted file (e.g.,
    default simulation configuration file).
    '''

    # ..................{ PROPERTIES                        }..................
    # Subclasses are required to implement these abstract class properties.

    @abstractclassproperty_readonly
    def conf_default_filename(cls) -> str:
        '''
        Absolute filename of the low-level **default YAML-formatted file**
        (i.e., file specific to this subclass, typically containing sane
        defaults of general interest to a wide audience).
        '''

        pass

    # ..................{ COPIERS                           }..................
    @type_check
    def copy_default(self, trg_conf_filename: str, **kwargs) -> None:
        '''
        Copy the default YAML-formatted file specific to this subclass to the
        YAML-formatted file with the passed filename *and* all external
        resources from the directory of the former into that of the latter.

        This method effectively implements the "New..." GUI metaphor.

        Parameters
        ----------
        trg_conf_filename : str
            Absolute or relative filename of the target YAML-formatted file to
            safely copy from.

        All remaining keyword arguments are passed to the :meth:`save` method.

        See Also
        ----------
        :meth:`copy`
            Further details.
        '''

        # Thus spake BETSEthustra.
        self.copy(
            src_conf_filename=self.conf_default_filename,
            trg_conf_filename=trg_conf_filename,
            **kwargs
        )
