#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
YAML-formatted simulation configuration file functionality.
'''

#FIXME: Refactor all functions defined below into methods of the
#"betse.science.config.wrapper.SimConfigWrapper" class.
#FIXME: Validate the versions of loaded configuration files.

# ....................{ IMPORTS                            }....................
from betse.lib.yaml import yamls
from betse.util.io.log import logs
from betse.util.path import files, paths
from betse.util.type.types import type_check, MappingType

# ....................{ LOADERS                            }....................
@type_check
def read(config_filename: str) -> MappingType:
    '''
    Deserialize the passed YAML-formatted simulation configuration file into a
    dictionary; then, validate and return this dictionary.

    Parameters
    ----------
    config_filename : str
        Absolute or relative path of the source YAML file to be deserialized.

    Returns
    ----------
    MappingType
        Dictionary deserialized from this file.
    '''

    # Load this dictionary from this YAML file.
    config = yamls.load(config_filename)

    #FIXME: Implement me *AFTER* the structure of such file settles down a tad.
    # Validate the contents of this file.

    # Return this dictionary.
    return config


#FIXME: Fix docstring and code duplicated from above.
@type_check
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

    # Load this dictionary from this YAML file.
    config = yamls.load(config_filename)

    #FIXME: Implement me *AFTER* the structure of such file settles down a tad.
    # Validate the contents of this file.

    # Return this dictionary.
    return config

# ....................{ WRITERS                            }....................
@type_check
def write(config_filename: str, config: MappingType) -> None:
    '''
    Serialize the passed dictionary to the passed YAML-formatted simulation
    configuration file.

    Parameters
    ----------
    config_filename : str
        Absolute or relative path of the target YAML file to be written.
    config : MappingType
        Dictionary to serialize to this file.

    Raises
    ----------
    BetseFileException
        If this file already exists.
    '''

    # Validate this file *BEFORE* writing this file.
    _write_check(config_filename)

    # Save this dictionary to this YAML file.
    yamls.save(config, config_filename)

# ....................{ WRITERS ~ default                  }....................
@type_check
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

    # Basename and filetype of this file.
    config_basename = paths.get_basename(config_filename)
    config_filetype = paths.get_filetype_undotted_or_none(config_basename)

    # If this file already exists, fail.
    files.die_if_file(config_filename)

    # If this filename is *NOT* suffixed by either ".yml" or ".yaml", log a
    # warning.
    if config_filetype not in ('yaml', 'yml'):
        logs.log_warning(
            'File "%s" filetype "%s" not "yaml" or "yml".',
            config_basename, config_filetype)
