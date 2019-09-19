#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level support facilities for Yet Another Markup Language (YAML), the
file format encapsulating most input and output data for this application.
'''

#FIXME: Consider all or most of the related "yamlrepr" submodule back to
#"ruamel.yaml" -- especially the non-trivial Numpy-to-YAML-type conversions.

# ....................{ IMPORTS                           }....................
from betse.util.io import iofiles
from betse.util.io.error.errwarning import ignoring_warnings
from betse.util.io.log import logs
from betse.util.path import pathnames
from betse.util.type.contexts import noop_context
from betse.util.type.obj import objects
from betse.util.type.types import (
    type_check, MappingOrSequenceTypes, StrOrNoneTypes)
from ruamel import yaml as ruamel_yaml
from ruamel.yaml.error import MantissaNoDotYAML1_1Warning

# ....................{ GLOBALS                           }....................
YAML_FILETYPES = {'yaml', 'yml',}
'''
Set of all YAML-compliant filetypes.
'''

# ....................{ LOADERS                           }....................
@type_check
def load(
    # Mandatory parameters.
    filename: str,

    # Optional parameters.
    yaml_version: StrOrNoneTypes = None,
) -> MappingOrSequenceTypes:
    '''
    Load (i.e., open and read, deserialize) and return the contents of the
    YAML-formatted file with the passed path as either a dictionary or list
    via the active YAML implementation.

    Parameters
    ----------
    filename : str
        Absolute or relative filename of the YAML-formatted file to be loaded.
    yaml_version: StrOrNoneTypes
        Version of the YAML specification (e.g., 1.1, 1.2) to forcefully assume
        this file to be compliant with, overriding any version directive
        prefacing this file (e.g., ``%YAML 1.2``). Defaults to ``None``, in
        which case the version directive prefacing this file is deferred to.

    Returns
    ----------
    MappingOrSequenceTypes
        Dictionary or list corresponding to the contents of this file.
    '''

    # If this filename has no YAML-compliant filetype, log a warning.
    _warn_unless_filetype_yaml(filename)

    # With this YAML file opened for character-oriented reading...
    with iofiles.reading_chars(filename) as yaml_file:
        # Safe roundtripping YAML parser.
        ruamel_parser = _make_ruamel_parser()

        # Context manager with which to load this file from this parser,
        # defaulting to a noop context manager.
        context_manager = noop_context()

        # If assuming this file to comply with a specific version of the
        # YAML specification...
        if yaml_version is not None:
            # Inform this parser of that.
            #
            # Sadly, this assignment is currently ignored by the subsequent
            # call to the load() method. See also this outstanding issue:
            #     https://bitbucket.org/ruamel/yaml/issues/320/yamlload-ignores-yamlversion
            ruamel_parser.version = yaml_version

            # If assuming this file to *NOT* comply with version 1.1 of this
            # specification, ignore YAML 1.1-specific warnings. Ideally, the
            # prior assignment would suffice; since it does *NOT*, we have no
            # recourse but to ignore these warnings explicitly.
            if yaml_version != '1.1':
                context_manager = ignoring_warnings(
                    MantissaNoDotYAML1_1Warning)

        # Load and return the contents of this file with this context manager.
        with context_manager:
            return ruamel_parser.load(yaml_file)

# ....................{ SAVERS                            }....................
@type_check
def save(
    # Mandatory parameters.
    container: MappingOrSequenceTypes,
    filename: str,

    # Optional parameters.
    is_overwritable: bool = False,
) -> None:
    '''
    Save (i.e., open and write, serialize) the passed dictionary or list to the
    YAML-formatted file with the passed path via the active YAML
    implementation.

    Parameters
    ----------
    container: MappingOrSequenceTypes
        Dictionary or list to be written as the contents of this file.
    filename : str
        Absolute or relative filename of the YAML-formatted file to be saved.
    is_overwritable : optional[bool]
        Either:

        * ``True`` if this function may silently overwrite this file when this
          file already exists.
        * ``False`` if this function should instead raise an exception when
          this file already exists.

        Defaults to ``False`` for safety.
    '''

    # If this filename has no YAML-compliant filetype, log a warning.
    _warn_unless_filetype_yaml(filename)

    # With this YAML file opened for character-oriented writing...
    with iofiles.writing_chars(
        filename=filename, is_overwritable=is_overwritable) as yaml_file:
        # Fully-qualified name of the module defining this container subclass.
        container_class_module_name = objects.get_class_module_name_qualified(
            container)

        # If this container is *NOT* a "ruamel.yaml"-specific object returned
        # by a prior call to the load() function, log a non-fatal warning.
        # While this edge case does *NOT* constitute a fatal warning, it does
        # disable roundtripped preservation of comments and whitespace
        # contained in the original YAML file -- the principal motivation for
        # preferring "ruamel.yaml" versus PyYAML and friends (e.g., oyaml).
        if not container_class_module_name.startswith('ruamel.'):
            logs.log_warning(
                'Non-"ruamel.yaml" type "%s.%s" not roundtrippable.',
                container_class_module_name, container.__class__.__name__)

        # Safe roundtripping YAML parser.
        ruamel_parser = _make_ruamel_parser()

        # Save this container to this YAML file.
        ruamel_parser.dump(container, yaml_file)

# ....................{ MAKERS                            }....................
def _make_ruamel_parser() -> ruamel_yaml.YAML:
    '''
    Safe roundtripping :mod:`ruamel.yaml` parser, where:

    * "Safe" implies this parser ignores all pragmas in this YAML file
      instructing parsers to construct arbitrary Python objects, whose YAML
      syntax is of the form: ``!!python/object:module.name { ... state ... }``.
    * "Roundtripping" implies the object deserialized from this YAML file
      preserves *all* comments and whitespace of this file, such that
      re-serializing this object back to disk preserves these substrings.
    '''

    # Avoid circular import dependencies.
    from betse.lib.yaml import yamlrepr

    # Safe roundtripping YAML parser.
    ruamel_parser = ruamel_yaml.YAML(
        # Type of YAML parser to produce. Note that, by design, the
        # "rt" (i.e., roundtripping) parser is *ALWAYS* guaranteed to be safe.
        typ='rt',
    )

    # Permit this parser to roundtrip Unicode characters.
    ruamel_parser.allow_unicode = True

    # Coerce this parser into serializing nested collections as blocks rather
    # than inline "flow": e.g.,
    #
    #     # ...like this.
    #     b:
    #       c: 3
    #       d: 4
    #
    #     # ...not this.
    #     b: {c: 3, d: 4}
    ruamel_parser.default_flow_style = False

    # Locally add all representers required by this application to this parser.
    yamlrepr.add_representers(ruamel_parser.representer)

    # Return this parser.
    return ruamel_parser

# ....................{ PRIVATE ~ warners                 }....................
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
