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
from betse import pathtree
from betse.lib.yaml import yamls
from betse.util.io.log import logs
from betse.util.path import dirs, files, pathnames
from betse.util.type.types import type_check, MappingType

# ....................{ EXCEPTIONS                         }....................
@type_check
def die_unless_writable(conf_filename: str) -> None:
    '''
    Raise an exception unless the passed path is suitable for writing a
    YAML-formatted simulation configuration file to.

    Specifically:

    * If this file already exists, an exception is raised for safety.
    * If this file's filetype is neither ``.yaml`` nor ``.yml``, a non-fatal
      warning is logged.

    Parameters
    ----------------------------
    conf_filename : str
        Absolute or relative path of the target YAML file to be validated.

    Raises
    ----------
    BetseFileException
        If this file already exists.
    '''

    # Basename and filetype of this file.
    config_basename = pathnames.get_basename(conf_filename)
    config_filetype = pathnames.get_filetype_undotted_or_none(config_basename)

    # If this file already exists, fail.
    files.die_if_file(conf_filename)

    # If this filename is suffixed by neither ".yml" nor ".yaml", log a warning.
    if config_filetype not in ('yaml', 'yml'):
        logs.log_warning(
            'Config file "%s" filetype "%s" neither "yaml" nor "yml".',
            config_basename, config_filetype)

# ....................{ LOADERS                            }....................
#FIXME: All functions defined in the "LOADERS" and "WRITERS" sections below are
#patently silly and, ideally, should be replaced by the following non-silly
#alternative:
#
#* Define a new "SimConfSerializableABC" superclass subclassing "SimConfABC" in
#  the "confabc" submodule.
#* In this class:
#  * Define a private _read() method identical to read() below.
#  * Define a private _write() method identical to write() below.
#* Subclass the "betse.science.parameters.Parameters" class from this new
#  "SimConfSerializableABC" superclass.
#* Likewise for all similar top-level classes currently calling these methods
#  such as "MasterOfMolecules", perhaps?

@type_check
def read(conf_filename: str) -> MappingType:
    '''
    Deserialize the passed YAML-formatted simulation configuration file into a
    dictionary; then, validate and return this dictionary.

    Parameters
    ----------
    conf_filename : str
        Absolute or relative path of the source YAML file to be deserialized.

    Returns
    ----------
    MappingType
        Dictionary deserialized from this file.
    '''

    # Load this dictionary from this YAML file.
    config = yamls.load(conf_filename)

    #FIXME: Implement me *AFTER* the structure of such file settles down a tad.
    # Validate the contents of this file.

    # Return this dictionary.
    return config


#FIXME: Fix docstring and code duplicated from above.
@type_check
def read_metabo(conf_filename: str) -> dict:
    '''
    Deserialize the passed YAML-formatted simulation configuration file into a
    dictionary; then, validate and return this dictionary.

    Parameters
    ----------
    conf_filename : str
        Absolute or relative path of the source YAML file to be deserialized.

    Returns
    ----------
    dict
        Dictionary deserialized from this file.
    '''

    # Load this dictionary from this YAML file.
    config = yamls.load(conf_filename)

    #FIXME: Implement me *AFTER* the structure of such file settles down a tad.
    # Validate the contents of this file.

    # Return this dictionary.
    return config

# ....................{ WRITERS                            }....................
@type_check
def write(conf_filename: str, config: MappingType) -> None:
    '''
    Serialize the passed dictionary to the passed YAML-formatted simulation
    configuration file.

    Parameters
    ----------
    conf_filename : str
        Absolute or relative path of the target YAML file to be written.
    config : MappingType
        Dictionary to serialize to this file.

    Raises
    ----------
    BetseFileException
        If this file already exists.
    '''

    # Validate this file *BEFORE* writing this file.
    die_unless_writable(conf_filename)

    # Save this dictionary to this YAML file.
    yamls.save(config, conf_filename)


@type_check
def write_default(conf_filename: str) -> None:
    '''
    Write a default YAML simulation configuration to the file with the passed
    path *and* recursively copy all external resources (e.g., images) required
    by this configuration into this file's parent directory.

    The resulting configuration will be usable as is with all high-level
    functionality requiring a valid configuration file (e.g., ``betse seed``).

    Parameters
    ----------
    conf_filename : str
        Absolute or relative path of the target YAML file to be written.

    Raises
    ----------
    BetseFileException
        If this file already exists.
    '''

    # Announce the ugly shape of things to come.
    logs.log_info('Writing default simulation configuration...')

    # Validate this file *BEFORE* writing this file.
    die_unless_writable(conf_filename)

    # Copy the default source simulation configuration file to this target file.
    files.copy(pathtree.get_sim_config_default_filename(), conf_filename)

    # Source directory containing the default simulation configuration.
    src_dirname = pathtree.get_data_yaml_dirname()

    # Note that the simple solution of recursively copying this source directory
    # into the parent directory of the passed target file (e.g., by calling
    # "dirs.copy(src_dirname, pathnames.get_dirname(conf_filename))") fails
    # for the following subtle reasons:
    #
    # * This target directory may be already exist, which dirs.copy() prohibits
    #   even if the directory is empty.
    # * This target configuration file basename may differ from that of the
    #   default source configuration file basename, necessitating a subsequent
    #   call to file.move().
    #
    # For each subdirectory of this directory...
    for src_subdirname in dirs.iter_subdirnames(src_dirname):
        # Recursively copy from this subdirectory into the target directory.
        dirs.copy_into_dir(
            src_dirname=src_subdirname,
            trg_dirname=pathnames.get_dirname(conf_filename),

            # Ignore all empty ".gitignore" files in all subdirectories of this
            # source directory. These files serve as placeholders instructing
            # Git to track their parent subdirectories, but otherwise serve no
            # purpose. Preserving these files only invites end user confusion.
            ignore_basename_globs=('.gitignore',),
        )
