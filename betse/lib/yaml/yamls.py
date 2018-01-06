#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level support facilities for Yet Another Markup Language (YAML), the
file format encapsulating most input and output data for this application.
'''

#FIXME: Consider contributing various portions of this submodule back to
#"ruamel.yaml" -- particularly the Numpy-type-to-YAML-native-type conversions.

# ....................{ IMPORTS                            }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To permit YAML implementations to be conditionally imported at
# application startup, no implementations (e.g., the top-level "yaml" package
# corresponding to PyYAML) are importable here.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

from betse import metadeps
from betse.exceptions import BetseYamlException
from betse.lib import libs
from betse.util.io import iofiles
from betse.util.io.log import logs
from betse.util.path import pathnames
from betse.util.type.obj import objects
from betse.util.type.types import (
    type_check, FileType, MappingOrSequenceTypes, ModuleType)
from enum import Enum

# ....................{ ENUMS                              }....................
YamlType = Enum('YamlType', ('PyYAML', 'ruamel_yaml',))
'''
Enumeration of all third-party YAML implementations supported by this submodule.
'''


YAML_TYPE_ACTIVE = None
'''
:data:`YamlType` member corresponding to the currently active third-party YAML
implementation if any *or* ``None`` if this submodule has yet to be initialized.
'''

# ....................{ GLOBALS                            }....................
FILETYPES = {'yaml', 'yml',}
'''
Set of all YAML-compliant filetypes.
'''

# ....................{ WARNERS                            }....................
@type_check
def warn_unless_filetype_yaml(filename: str) -> None:
    '''
    Log a non-fatal warning unless the passed filename has a YAML-compliant
    filetype (i.e., is suffixed by either ``.yaml`` or ``.yml``).

    Parameters
    ----------
    filename : str
        Absolute or relative path of this file.
    '''

    # Filetype of this file.
    filetype = pathnames.get_filetype_undotted_or_none(filename)

    # If this filetype is *NOT* YAML-compliant...
    if filetype not in FILETYPES:
        # If this file has *NO* filetype, log an appropriate warning.
        if filetype is None:
            logs.log_warning('YAML file "%s" has no filetype.', filename)
        # Else, this file has a filetype. Log an appropriate warning.
        else:
            logs.log_warning(
                'YAML file "%s" filetype "%s" neither "yaml" nor "yml".',
                filename, filetype)

# ....................{ LOADERS                            }....................
@type_check
def load(filename: str) -> MappingOrSequenceTypes:
    '''
    Load (i.e., open and read, deserialize) and return the contents of the
    YAML-formatted file with the passed path as either a dictionary or list
    via the active YAML implementation.

    Parameters
    ----------
    filename : str
        Absolute or relative path of this file.

    Returns
    ----------
    MappingOrSequenceTypes
        Dictionary or list corresponding to the contents of this file.
    '''

    # If this filename has no YAML-compliant filetype, log a warning.
    warn_unless_filetype_yaml(filename)

    # With this YAML file opened for character-oriented reading...
    with iofiles.reading_chars(filename) as yaml_file:
        # Load this YAML file via the active YAML implementation.
        if YAML_TYPE_ACTIVE is YamlType.PyYAML:
            return _load_pyyaml(yaml_file)
        elif YAML_TYPE_ACTIVE is YamlType.ruamel_yaml:
            return _load_ruamel(yaml_file)
        else:
            raise BetseYamlException(
                'YAML type "{}" unrecognized.'.format(YAML_TYPE_ACTIVE))


@type_check
def _load_pyyaml(yaml_file: FileType) -> MappingOrSequenceTypes:
    '''
    Load and return the contents of the YAML-formatted file with the passed
    readable file handle via the PyYAML parser (i.e., the :mod:`yaml` package).
    '''

    # PyYAML, validated and imported dynamically for safety.
    pyyaml = _import_pyyaml()

    # Load and return the contents of this YAML file.
    return pyyaml.load(yaml_file)


@type_check
def _load_ruamel(yaml_file: FileType) -> MappingOrSequenceTypes:
    '''
    Load and return the contents of the YAML-formatted file with the passed
    readable file handle via the :mod:`ruamel.yaml` parser.
    '''

    # Safe roundtripping YAML parser.
    ruamel_parser = _make_ruamel_parser()

    # Load and return the contents of this YAML file.
    return ruamel_parser.load(yaml_file)

# ....................{ SAVERS                             }....................
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
    YAML-formatted file with the passed path via the active YAML implementation.

    Parameters
    ----------
    container: MappingOrSequenceTypes
        Dictionary or list to be written as the contents of this file.
    filename : str
        Absolute or relative path of this file.
    is_overwritable : optional[bool]
        ``True`` if overwriting this file when this file already exists *or*
        ``False`` if raising an exception when this file already exists.
        Defaults to ``False`` for safety.
    '''

    # If this filename has no YAML-compliant filetype, log a warning.
    warn_unless_filetype_yaml(filename)

    # With this YAML file opened for character-oriented writing...
    with iofiles.writing_chars(
        filename=filename, is_overwritable=is_overwritable) as yaml_file:
        # Save this YAML file via the active YAML implementation.
        if YAML_TYPE_ACTIVE is YamlType.PyYAML:
            return _save_pyyaml(container, yaml_file)
        elif YAML_TYPE_ACTIVE is YamlType.ruamel_yaml:
            return _save_ruamel(container, yaml_file)
        else:
            raise BetseYamlException(
                'YAML type "{}" unrecognized.'.format(YAML_TYPE_ACTIVE))

@type_check
def _save_pyyaml(
    container: MappingOrSequenceTypes, yaml_file: FileType) -> None:
    '''
    Save the passed container to the YAML-formatted file with the passed
    writable file handle via the PyYAML parser (i.e., the :mod:`yaml` package).
    '''

    # PyYAML, validated and imported dynamically for safety.
    pyyaml = _import_pyyaml()

    # Save this container to this YAML file.
    pyyaml.dump(
        data=container,
        stream=yaml_file,

        # Permit this parser to roundtrip Unicode characters.
        allow_unicode=True,
        encoding=None,

        # Coerce this parser into serializing nested collections as blocks
        # rather than inline "flow". See _make_ruamel_parser() for details.
        default_flow_style=False,
    )


@type_check
def _save_ruamel(
    container: MappingOrSequenceTypes, yaml_file: FileType) -> None:
    '''
    Save the passed container to the YAML-formatted file with the passed
    writable file handle via the :mod:`ruamel.yaml` parser.
    '''

    # Fully-qualified name of the module defining this container's subclass.
    container_class_module_name = objects.get_class_module_name(container)

    # If this container is *NOT* a "ruamel.yaml"-specific object returned by a
    # prior call to the _load_ruamel() function, log a non-fatal warning. While
    # this edge case does *NOT* constitute a fatal warning, it does disable
    # roundtripped preservation of comments and whitespace contained in the
    # original YAML file -- the principal motivation for preferring
    # "ruamel.yaml" over PyYAML.
    if not container_class_module_name.startswith('ruamel.'):
        logs.log_warning(
            'Non-"ruamel.yaml" type "%s.%s" not roundtrippable.',
            container_class_module_name, container.__class__.__name__)

    # Safe roundtripping YAML parser.
    ruamel_parser = _make_ruamel_parser()

    # Save this container to this YAML file.
    ruamel_parser.dump(container, yaml_file)

# ....................{ INITIALIZERS                       }....................
def init() -> None:
    '''
    Initialize the currently active YAML implementation (e.g.,
    :mod:`ruamel.yaml`, PyYAML).

    This function selects the first such implementation importable under the
    active Python interpreter from the following list (in descending order of
    preference): :mod:`ruamel.yaml`, PyYAML.
    '''

    # Globals assigned to below.
    global YAML_TYPE_ACTIVE

    # Log the active YAML implementation.
    logs.log_debug(
        'Initializing preferred YAML binding "%s"...',
        metadeps.RUNTIME_MANDATORY_YAML_PROJECT_NAME)

    # Convert the non-typesafe setuptools-specific project name of this
    # implementation into a typesafe enumeration member.
    if metadeps.RUNTIME_MANDATORY_YAML_PROJECT_NAME == 'PyYAML':
        YAML_TYPE_ACTIVE = YamlType.PyYAML
    elif metadeps.RUNTIME_MANDATORY_YAML_PROJECT_NAME == 'ruamel.yaml':
        YAML_TYPE_ACTIVE = YamlType.ruamel_yaml
    else:
        raise BetseYamlException(
            'YAML project name "{}" unrecognized.'.format(
                metadeps.RUNTIME_MANDATORY_YAML_PROJECT_NAME))

    # If the active YAML implementation is PyYAML, initialize PyYAML.
    #
    # Since "ruamel.yaml" prefers a modern object-oriented API and hence is
    # locally initialized on object construction rather than globally on module
    # importation, note that "ruamel.yaml" requires no such initialization.
    if YAML_TYPE_ACTIVE is YamlType.PyYAML:
        _init_pyyaml()


def _init_pyyaml() -> None:
    '''
    Initialize PyYAML as the currently active YAML implementation.
    '''

    # Avoid circular import dependencies.
    from betse.lib.yaml import yamlrepr

    # PyYAML, validated and imported dynamically for safety.
    pyyaml = _import_pyyaml()

    # Globally add all representers required by this application to PyYAML.
    yamlrepr.add_representers(pyyaml)

    # Monkeypatch the following PyYAML method to eliminate Numpy-specific
    # warnings. Note that "ruamel.yaml" has since addressed all
    # Numpy-specific warnings and thus does *NOT* require this monkeypatch.
    pyyaml.representer.SafeRepresenter.ignore_aliases = _pyyaml_ignore_aliases

# ....................{ IMPORTERS                          }....................
def _import_pyyaml() -> ModuleType:
    '''
    Dynamically validate, import, and return the top-level module for PyYAML.
    '''

    return libs.import_runtime_optional('PyYAML')


def _import_ruamel() -> ModuleType:
    '''
    Dynamically validate, import, and return the top-level module for
    :mod:`ruamel.yaml`.
    '''

    return libs.import_runtime_optional('ruamel.yaml')

# ....................{ MAKERS                             }....................
def _make_ruamel_parser() -> 'ruamel.yaml.YAML':
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

    # "ruamel.yaml", validated and imported dynamically for safety.
    ruamel_yaml = _import_ruamel()

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

# ....................{ MONKEYPATCHES                      }....................
def _pyyaml_ignore_aliases(self, data) -> bool:
    '''
    PyYAML-specific :meth:`yaml.representer.SafeRepresenter.ignore_aliases`
    method monkeypatched to eliminate Numpy-specific warnings.

    This method has been refactored to eliminate future warnings resembling:

        /usr/lib64/python3.4/site-packages/yaml/representer.py:135:
        FutureWarning: comparison to `None` will result in an elementwise object
        comparison in the future.
          if data in [None, ()]:
    '''

    # Refactor the default test of "if data in [None, ()]:" in the stock
    # ignore_aliases() method to avoid the future warning detailed above. To do
    # so safely, note that testing the passed "data" parameter against:
    #
    # * The singleton "None" object with the "is" operator is required.
    # * The non-singleton () object with the "==" rather than "is" operator is
    #   required. Bizarrely, despite the empty tuple being frozen, Python
    #   documentation notes that "...two occurrences of the empty tuple may or
    #   may not yield the same object". Hence, the comparison "data is ()" could
    #   unsafely return "False" even when "data" is the empty tuple!
    if data is None or data == ():
        return True
    if isinstance(data, (str, bytes, bool, int, float)):
        return True
