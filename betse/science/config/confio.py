#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Simulation configuration file input and output (I/O) facilities.
'''

#FIXME: Validate the versions of loaded configuration files.

# ....................{ IMPORTS                            }....................
from betse import pathtree
from betse.lib.yaml import yamls
from betse.util.io.log import logs
from betse.util.path import dirs, files, pathnames
from betse.util.type.types import type_check  #, GeneratorType

# ....................{ READERS                            }....................
#FIXME: This utility function is patently silly and should be replaced by the
#following non-silly alternative:
#
#* Subclass all top-level classes currently calling this function (e.g.,
#  "MasterOfMolecules") from the the "betse.lib.yaml.yamlabc.YamlFileABC" class.
#* Replace all calls to this function with calls to "self.read(conf_filename)".

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
def write_default(conf_filename: str, is_overwritable: bool = False) -> None:
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
    is_overwritable : optional[bool]
        For any path to be written by this function that already exists when
        this boolean is:
        * ``True``, that path is silently overwritten.
        * ``False``, an exception is raised.
        Defaults to ``False``.

    Raises
    ----------
    BetseFileException
        If this file already exists *and* ``is_overwritable`` is ``False``.
    '''

    # Announce the ugly shape of things to come.
    logs.log_info('Writing default simulation configuration...')

    # Copy the default source simulation configuration file to this target file.
    files.copy(
        src_filename=pathtree.get_sim_config_default_filename(),
        trg_filename=conf_filename,
        is_overwritable=is_overwritable,
    )

    # Source directory containing the default simulation configuration.
    src_dirname = pathtree.get_data_yaml_dirname()

    # Target directory to be copied into, defined to be either the parent
    # directory of the passed path if this path has a dirname or the current
    # working directory otherwise.
    trg_dirname = pathnames.get_dirname_or_cwd(conf_filename)

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
            trg_dirname=trg_dirname,
            is_overwritable=is_overwritable,

            # Ignore all empty ".gitignore" files in all subdirectories of this
            # source directory. These files serve as placeholders instructing
            # Git to track their parent subdirectories, but otherwise serve no
            # purpose. Preserving these files only invites end user confusion.
            ignore_basename_globs=('.gitignore',),
        )
