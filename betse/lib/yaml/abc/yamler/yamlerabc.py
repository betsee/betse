#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Abstract base classes of all YAML-backed simulation configuration subclasses as
well as functionality pertaining to such classes.
'''

# ....................{ IMPORTS                           }....................
from abc import ABCMeta, abstractmethod
# from betse.exceptions import BetseYamlException
from betse.lib import libs
from betse.util.io import iofiles
from betse.util.io.log import logs
from betse.util.path import pathnames
from betse.util.type.obj import objects
from betse.util.type.types import (
    type_check, FileType, MappingOrSequenceTypes,)

# ....................{ SUPERCLASSES                      }....................
class YamlerABC(object, metaclass=ABCMeta):
    '''
    Abstract base class of all **YAML implementation wrapper** (i.e.,
    portable high-level object wrapping a low-level third-party YAML
    (de)serialization library) subclasses.

    This superclass enables callers to transparently switch between third-party
    YAML implementations (e.g., PyYAML, :mod:`oyaml`, :mod:`ruamel.yaml`),
    whose public APIs are increasingly incompatible with one another and hence
    *not* safely interchangeable as is.

    Attributes
    ----------
    '''

    # ..................{ INITIALIZERS                      }..................
    # def __init__(self) -> None:
    #     '''
    #     Initialize this YAML implementation wrapper.
    #     '''
    #
    #     pass

    # ..................{ SUBCLASS                          }..................
    # Subclasses are required to implement the following abstract methods.

    @abstractmethod
    def _load_opened(self, yaml_file: FileType) -> MappingOrSequenceTypes:
        '''
        Load (i.e., open and read, deserialize) and return the contents of the
        YAML-formatted file with the passed readable open file handle as either
        a dictionary or list via this YAML implementation.

        Parameters
        ----------
        yaml_file : FileType
            Readable open file handle of the YAML file to be loaded.

        Returns
        ----------
        MappingOrSequenceTypes
            Dictionary or list corresponding to the contents of this file.

        See Also
        ----------
        :meth:`load`
            Higher-level implementation-agnostic method wrapping this
            lower-level implementation-specific method.
        '''

        pass

    # ..................{ LOADERS                           }..................
    @type_check
    def load(self, filename: str) -> MappingOrSequenceTypes:
        '''
        Load (i.e., open and read, deserialize) and return the contents of the
        YAML-formatted file with the passed filename as either a dictionary or
        list with this YAML implementation.

        Parameters
        ----------
        filename : str
            Absolute or relative filename of the YAML file to be loaded.

        Returns
        ----------
        MappingOrSequenceTypes
            Dictionary or list corresponding to the contents of this file.
        '''

        # If this filename has no YAML-compliant filetype, log a warning.
        _warn_unless_filetype_yaml(filename)

        # With this YAML file temporarily opened for character-oriented reads,
        # load this YAML file via the active YAML implementation.
        with iofiles.reading_chars(filename) as yaml_file:
            return self._load_opened(yaml_file)

# ....................{ SUBCLASSES                        }....................
class YamlerPyYAML(YamlerABC):
    '''
    **PyYAML implementation wrapper** (i.e.,
    portable high-level object wrapping the low-level third-party :mod:`yaml`
    (de)serialization library).
    '''

    # ..................{ SUPERCLASS                        }..................
    @type_check
    def _load_opened(self, yaml_file: FileType) -> MappingOrSequenceTypes:

        # PyYAML, validated and imported dynamically for safety.
        pyyaml = libs.import_runtime_optional('PyYAML')

        # PyYAML function loading and returning the contents of this YAML file,
        # dynamically obtained in a version-agnostic manner preserving backward
        # compatibility with obsolete (albeit commonplace) versions of PyYAML.
        # Specifically:
        #
        # * If this version of PyYAML defines the full_load() function and is
        #   thus sufficiently modern, this function is strongly preferred. This
        #   function implements the default PyYAML loader (i.e.,
        #   "yaml.FullLoader"), the only PyYAML loader both supporting the full
        #   YAML specification *AND* avoiding arbitrary code execution.
        #
        #   Note this is effectively syntactic sugar for the following call:
        #       return pyyaml.load(input, Loader=pyyaml.FullLoader)
        #
        #   Note that the pyyaml.load() function defaults to this loader but
        #   emits a non-fatal (albeit frightening) deprecation warning on doing
        #   so. Ergo, this loader *MUST* effectively be explicitly specified.
        # * Else, this version of PyYAML is obsolete. In this case, the load()
        #   function known to suffer security vulnerabilities is deferred to.
        #
        # For exhaustive details, see also:
        #     https://github.com/yaml/pyyaml/wiki/PyYAML-yaml.load(input)-Deprecation
        pyyaml_loader = (
            objects.get_callable_or_none(
                obj=pyyaml, callable_name='full_load') or
            objects.get_callable_or_none(
                obj=pyyaml, callable_name='load'))

        # Return the contents of this YAML file with this loader.
        return pyyaml_loader(yaml_file)

# ....................{ WARNERS                           }....................
@type_check
def _warn_unless_filetype_yaml(filename: str) -> None:
    '''
    Log a non-fatal warning unless the passed filename has a YAML-compliant
    filetype (i.e., is suffixed by either ``.yaml`` or ``.yml``).

    Parameters
    ----------
    filename : str
        Absolute or relative filename of this file.
    '''

    # Avoid circular import dependencies.
    from betse.lib.yaml.yamls import YAML_FILETYPES

    # Filetype of this file.
    filetype = pathnames.get_filetype_undotted_or_none(filename)

    # If this filetype is *NOT* YAML-compliant...
    if filetype not in YAML_FILETYPES:
        # If this file has *NO* filetype, log an appropriate warning.
        if filetype is None:
            logs.log_warning('YAML file "%s" has no filetype.', filename)
        # Else, this file has a filetype. Log an appropriate warning.
        else:
            logs.log_warning(
                'YAML file "%s" filetype "%s" neither "yaml" nor "yml".',
                filename, filetype)
