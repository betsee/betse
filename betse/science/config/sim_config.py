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
from betse import pathtree
from betse.util.io.log import logs
from betse.util.path import dirs, files, paths
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
    with files.open_for_text_reading(config_filename) as yaml_file:
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
    with files.open_for_text_reading(config_filename) as yaml_file:
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
    BetseExceptionFile
        If this file already exists.
    '''
    assert types.is_str_nonempty(config_filename), (
        types.assert_not_str_nonempty(config_filename, 'Filename'))
    assert types.is_mapping(config), types.assert_not_mapping(config)

    # Validate this file *BEFORE* writing this file.
    _write_check(config_filename)

    # Open this filename for writing and...
    with files.open_for_text_writing(config_filename) as yaml_file:
        # String serialized from this dictionary.
        config_dump = yaml.dump(
            config,
            allow_unicode=True, default_flow_style=False, encoding=None)

        # Write this string to this file.
        yaml_file.write(config_dump)

# ....................{ WRITERS ~ default                  }....................
#FIXME: Refactor this function as follows (in order):
#
#* Create a new "betse.science.config.default" submodule.
#* Shift this and all related private functions of this submodule into that
#  submodule.
#* Rename this function to merely write() and likewise for all other functions
#  moved into that submodule.
#* Refactor this function to internally create and return a "SimConfigWrapper"
#  instance. Since this function is only called once elsewhere in the codebase,
#  this should be relatively simple. Do so now, however, before matters become
#  more... entangled.
def write_default(config_filename: str) -> None:
    '''
    Write a default YAML simulation configuration to the file with the passed
    path _and_ recursively copy all external resources (e.g., images) required
    by this configuration into this file's parent directory.

    The resulting configuration will be usable as is with all high-level `betse`
    functionality requiring a valid configuration file (e.g., `betse world`).

    Modifications
    ----------
    For usability, this method modifies the contents of the written (but _not_
    original) file as follows:

    * The `plot after solving` option in the `results options` section is set to
      `False`. Ideally, this prevents hapless end-users from drowning under an
      intimidating deluge of plot windows irrelevant to general usage.

    Parameters
    ----------
    config_filename : str
        Absolute or relative path of the target YAML file to be written.

    Raises
    ----------
    BetseExceptionFile
        If this file already exists.
    '''

    # Announce the ugly shape of things to come.
    logs.log_info('Writing default simulation configuration.')

    # Validate this file *BEFORE* creating this file's parent directory if
    # needed *BEFORE* creating this file.
    _write_check(config_filename)
    _write_default_dir(config_filename)
    _write_default_file(config_filename)


def _write_default_dir(config_filename: str) -> None:
    '''
    Recursively copy all **external resources** (i.e., data-specific
    subdirectories and files) referenced by the default YAML configuration file
    into the parent directory of the passed path.

    If not currently found, this directory will also be created.

    Parameters
    ----------
    config_filename : str
        Absolute or relative path of the target YAML file to be written. If this
        file has no dirname and hence is a pure basename (e.g.,
        `sim_config.yaml`), the parent directory to which resources are copied
        will be the current working directory (CWD).
    '''
    assert types.is_str_nonempty(config_filename), (
        types.assert_not_str_nonempty(config_filename, 'Filename'))

    # Parent directory of this file if any or the current directory otherwise.
    target_dirname = paths.get_dirname_or_current_dirname(config_filename)

    # Create this directory if needed.
    dirs.make_unless_dir(target_dirname)

    # For the absolute path of each source subdirectory containing assets
    # required by this file...
    for source_asset_dirname in pathtree.DATA_DEFAULT_ASSET_DIRNAMES:
        # Absolute path of the corresponding target subdirectory.
        target_asset_dirname = paths.join(
            target_dirname, paths.get_basename(source_asset_dirname))

        # If this subdirectory already exists, log a non-fatal warning.
        if dirs.is_dir(target_asset_dirname):
            logs.log_warning(
                'Ignoring existing subdirectory "{}".'.format(
                    target_asset_dirname))
        # Else, recursively copy the entirety of this source subdirectory into
        # this target subdirectory.
        else:
            dirs.copy(source_asset_dirname, target_asset_dirname)


def _write_default_file(config_filename: str) -> None:
    '''
    Write a default YAML simulation configuration file to the passed path,
    assumed to _not_ already exist.

    Parameters
    ----------
    config_filename : str
        Absolute or relative path of the YAML file to be written.
    '''
    assert types.is_str_nonempty(config_filename), (
        types.assert_not_str_nonempty(config_filename, 'Filename'))

    #FIXME: Ideally, we should be using ruamel.yaml to munge YAML data in a
    #well-structured and hence sane manner rather than the admittedly crude
    #"sed"-like approach leveraged below. Unfortunately, given the complex
    #nature of this data, it's unclear whether or not ruamel.yaml would
    #adequately preserve the entirety of this data in a roundtrip manner. For
    #now, the "sed"-like approach prevails.

    # Write the default configuration to this file, modifying the latter with
    # "sed"-like global string substitution as detailed above.
    files.substitute_substrings(
        filename_source=pathtree.CONFIG_DEFAULT_FILENAME,
        filename_target=config_filename,
        substitutions=(
            # Prevent static plots from being displayed by default.
            (r'^(\s*plot after solving:\s+)True\b(.*)$', r'\1False\2'),
        ),
    )

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
    BetseExceptionFile
        If this file already exists.
    '''
    assert types.is_str_nonempty(config_filename), (
        types.assert_not_str_nonempty(config_filename, 'Filename'))

    # Basename and filetype of this file.
    config_basename = paths.get_basename(config_filename)
    config_filetype = paths.get_filetype(config_basename)

    # If this file already exists, fail.
    files.die_if_file(config_filename)

    # If this filename is *NOT* suffixed by either ".yml" or ".yaml", log a
    # warning.
    if not (config_filetype == 'yaml' or config_filetype == 'yml'):
        logs.log_warning(
            'File "{}" filetype "{}" not "yaml" or "yml".'.format(
                config_basename, config_filetype))
