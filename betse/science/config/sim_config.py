#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Simulation configuration in YAML format.
'''

#FIXME: Refactor all functions defined below into methods of the
#"betse.science.config.wrapper.SimConfigWrapper" class.
#FIXME: Validate the versions of loaded configuration files.

# ....................{ IMPORTS                            }....................
import yaml

from betse.util.io.log import logs
from betse.util.path import files, paths
from betse.util.type import types


# ....................{ LOADERS                            }....................
def read(config_filename: str) -> dict:
    '''
    Deserialize the passed YAML-formatted simulation configuration file into a
    dictionary; then, validate and return this dictionary.

    Parameters
    ----------
    config_filename : str
        Absolute or relative path of the source YAML file to be deserialized.

    Returns
    ----------
    dict
        Dictionary deserialized from this file.
    '''
    assert types.is_str_nonempty(config_filename), (
        types.assert_not_str_nonempty(config_filename, 'Filename'))

    # Dictionary deserialized from this file.
    config = None

    # Open this filename for reading and read this file into a dictionary.
    with files.read_chars(config_filename) as yaml_file:
        config = yaml.load(yaml_file)

    #FIXME: Implement me *AFTER* the structure of such file settles down a tad.
    # Validate the contents of this file.

    # Return this dictionary.
    return config


#FIXME: Fix docstring and code duplicated from above.
def read_metabo(config_filename: str) -> dict:
    '''
    Deserialize the passed YAML-formatted simulation configuration file into a
    dictionary; then, validate and return this dictionary.

    Parameters
    ----------
    config_filename : str
        Absolute or relative path of the source YAML file to be deserialized.

    Returns
    ----------
    dict
        Dictionary deserialized from this file.
    '''
    assert types.is_str_nonempty(config_filename), (
        types.assert_not_str_nonempty(config_filename, 'Filename'))

    # Dictionary deserialized from this file.
    config = None

    # Open this filename for reading and read this file into a dictionary.
    with files.read_chars(config_filename) as yaml_file:
        config = yaml.load(yaml_file)

    #FIXME: Implement me *AFTER* the structure of such file settles down a tad.
    # Validate the contents of this file.

    # Return this dictionary.
    return config

# ....................{ WRITERS                            }....................
def write(config_filename: str, config: dict) -> None:
    '''
    Serialize the passed dictionary to the passed YAML-formatted simulation
    configuration file.

    Parameters
    ----------
    config_filename : str
        Absolute or relative path of the target YAML file to be written.
    config : dict
        Dictionary to serialize to this file.

    Raises
    ----------
    BetseFileException
        If this file already exists.
    '''
    assert types.is_str_nonempty(config_filename), (
        types.assert_not_str_nonempty(config_filename, 'Filename'))
    assert types.is_mapping(config), types.assert_not_mapping(config)

    # Validate this file *BEFORE* writing this file.
    _write_check(config_filename)

    # Open this filename for writing and...
    with files.write_chars(config_filename) as yaml_file:
        # String serialized from this dictionary.
        config_dump = yaml.dump(
            config,
            allow_unicode=True, default_flow_style=False, encoding=None)

        # Write this string to this file.
        yaml_file.write(config_dump)

# ....................{ WRITERS ~ default                  }....................
def _write_check(config_filename: str) -> None:
    '''
    Validate the suitability of the passed path for use as a target YAML
    simulation configuration file.

    Specifically:

    * If this file already exists, an exception is raised for safety.
    * If this file's filetype is _not_ `.yaml`, a non-fatal warning is logged.

    Parameters
    ----------------------------
    config_filename : str
        Absolute or relative path of the target YAML file to be validated.

    Raises
    ----------
    BetseFileException
        If this file already exists.
    '''
    assert types.is_str_nonempty(config_filename), (
        types.assert_not_str_nonempty(config_filename, 'Filename'))

    # Basename and filetype of this file.
    config_basename = paths.get_basename(config_filename)
    config_filetype = paths.get_filetype_undotted_or_none(config_basename)

    # If this file already exists, fail.
    files.die_if_file(config_filename)

    # If this filename is *NOT* suffixed by either ".yml" or ".yaml", log a
    # warning.
    if not (config_filetype == 'yaml' or config_filetype == 'yml'):
        logs.log_warning(
            'File "{}" filetype "{}" not "yaml" or "yml".'.format(
                config_basename, config_filetype))
