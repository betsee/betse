#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Simulation configuration in YAML format.
'''

#FIXME: Validate the versions of loaded configuration files.

# ....................{ IMPORTS                            }....................
import yaml
from betse import pathtree
from betse.util.io import loggers
from betse.util.path import dirs, files, paths
from betse.util.type import types

# ....................{ LOADERS                            }....................
def load(config_filename: str) -> dict:
    '''
    Load and validate the passed YAML simulation configuration file.

    Parameters
    ----------------------------
    config_filename : str
        Absolute or relative path of the YAML file to be loaded.
    '''
    assert types.is_str_nonempty(config_filename), (
        types.assert_not_str_nonempty(config_filename, 'Filename'))

    # Contents of this file as a dictionary.
    config = None

    # Open the passed file for reading and read this file into a dictionary.
    with files.open_for_text_reading(config_filename) as yaml_file:
        config = yaml.load(yaml_file)

    #FIXME: Implement me *AFTER* the structure of such file settles down a tad.
    # Validate the contents of such file.

    # Return such dictionary.
    return config

# ....................{ WRITERS                            }....................
def write_default(config_filename: str) -> None:
    '''
    Write a default YAML simulation configuration file to the passed path _and_
    copy all external resources (e.g., images) referenced by this configuration
    into the parent directory of this path.

    The resulting configuration will be usable as is with all high-level `betse`
    functionality requiring a valid configuration file (e.g., `betse world`).

    Changes
    ----------------------------
    For usability, the contents of the written (but _not_ original)
    configuration file will be modified as follows:

    * The `plot after solving` option in the `results options` section will be
      forcefully set to `False`. Ideally, this prevents the hapless end user
      from drowning under an intimidating deluge of static plot windows
      irrelevant to general-purpose usage.

    Parameters
    ----------------------------
    config_filename : str
        Absolute or relative path of the target YAML file to be written.
    '''
    # Validate this file *BEFORE* creating this file's parent directory if
    # needed *BEFORE* creating this file.
    _write_default_check(config_filename)
    _write_default_dir(config_filename)
    _write_default_file(config_filename)


def _write_default_check(config_filename: str) -> None:
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
    '''
    assert types.is_str_nonempty(config_filename), (
        types.assert_not_str_nonempty(config_filename, 'Filename'))

    # Basename and filetype of this file.
    config_basename = paths.get_basename(config_filename)
    config_filetype = paths.get_filetype(config_basename)

    # Announce the ugly shape of things to come.
    loggers.log_info('Writing default simulation configuration.')

    # If this file already exists, fail. (For safety, we avoid silently
    # overwriting existing files.)
    files.die_if_file(config_filename)

    # If this filename is *NOT* suffixed by ".yaml", log a warning.
    if config_filetype != 'yaml':
        loggers.log_warning(
            'File "{}" filetype "{}" not "yaml".'.format(
                config_basename, config_filetype))


def _write_default_dir(config_filename: str) -> None:
    '''
    Recursively copy all **external resources** (i.e., data-specific
    subdirectories and files) referenced by the default YAML configuration file
    into the parent directory of the passed path.

    If not currently found, this directory will also be created.

    Parameters
    ----------------------------
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

    # Copy the source channel library to this directory.
    files.copy(pathtree.DATA_CHANNELS_YAML_FILENAME, target_dirname)

    # Absolute path of the target geometry subdirectory to be created below.
    target_geometry_dirname = paths.join(
        target_dirname, paths.get_basename(pathtree.DATA_GEOMETRY_DIRNAME))

    # If this directory already exists, log a non-fatal warning.
    if dirs.is_dir(target_geometry_dirname):
        loggers.log_warning(
            'Ignoring existing subdirectory "{}".'.format(
                target_geometry_dirname))
    # Else, copy the source geometry subdirectory to this directory.
    else:
        dirs.copy(pathtree.DATA_GEOMETRY_DIRNAME, target_geometry_dirname)


def _write_default_file(config_filename: str) -> None:
    '''
    Write a default YAML simulation configuration file to the passed path,
    assumed to _not_ already exist.

    Parameters
    ----------------------------
    config_filename : str
        Absolute or relative path of the YAML file to be written.
    '''
    assert types.is_str_nonempty(config_filename), (
        types.assert_not_str_nonempty(config_filename, 'Filename'))

    #FIXME: Ideally, we should be using ruamel.yaml to munge YAML data in a
    #well-structured and hence sane manner rather than the admittedly crude
    #"sed"-like approach leveraged below. Unfortunately, given the complex
    #nature of that data, it's unclear whether or not ruamel.yaml would
    #adequately preserve the entirety of that data in a roundtrip manner. For
    #now, the "sed"-like approach prevails.

    # Write the default configuration to this file, modifying the latter with
    # "sed"-like global string substitution as detailed above.
    files.substitute_substrings(
        filename_source = pathtree.CONFIG_DEFAULT_FILENAME,
        filename_target = config_filename,
        substitutions = (
            # Prevent static plots from being displayed by default.
            (r'^(\s*plot after solving:\s+)True\b(.*)$', r'\1False\2'),
        ),
    )

# --------------------( WASTELANDS                         )--------------------
    # Basename and filetype of this file.
    # config_basename = paths.get_basename(config_filename)
    # config_filetype = paths.get_filetype(config_basename)

    # target_dirname : str
    #     Absolute path of the target directory to which resources will be copied.
    # After loading this file, this function validates the contents of such file.

#Copy all external files referenced by such file to its parent directory.
    # Specifically, this function raises an exception if:
    #
    # *

    # # If such filename contains a dirname and hence is *NOT* a pure basename,
    # # create such file's parent directory if needed.
    # if dirname:
    #     dirs.make_unless_dir(dirname)
    # # Else, default such dirname to the current directory.
    # else:

    # loggers.log_info(
    #     'Copying file "%s" to "%s".',
    #     pathtree.CONFIG_DEFAULT_FILENAME, filename)
            # (r'(\s*turn all plots off)', r'\1'),
    # loggers.log_info(
    #     'Writing default simulation configuration to "{}".'.format(basename))

# ....................{ GETTERS                            }....................
# def get() -> str:
#     '''
#     Get the default YAML string for configuring simulations.
#
#     Such string is intended to be subsequently serialized to disk (e.g., as a
#     file suffixed by `.yaml`).
#     '''
#     return _TEMPLATE.format(
#         program_name = metadata.NAME,
#     )

# ....................{ CONSTANTS ~ private                }....................
# _TEMPLATE = '''
# %YAML 1.1
# ---
# # Default tissue simulation configuration.
# #
# # You are welcome to change any of the following settings.
# #
# # Pathnames
# # ----------
# # For portability, files and directories are configured below as relative rather
# # than absolute pathnames (e.g., as "sim_init.betse" rather than
# # "/Users/iamawesome/my_sim/sim_init.betse"). Such pathnames are relative to the
# # directory containing this file, allowing such directory to be safely moved
# # without breaking existing tissue simulations.
# #
# # YAML
# # ----------
# # As the "%YAML"-prefixed line above implies, this file is formatted according
# # to the YAML ("[Y]AML [A]in't [M]arkup [L]anguage") standard. YAML is a human-
# # readable data serialization format useful for configuration files -- like this
# # one. For details, see https://en.wikipedia.org/wiki/YAML.
#
# # -----------------------------------------------------------------------------
# # INITIALIZATION
# # -----------------------------------------------------------------------------
# # Configure simulation initialization.
# init:
#
#   # Files output by such initialization.
#   file:
#
#     # File caching the results of such initialization.
#     cache: sim_init.betse
#
# # -----------------------------------------------------------------------------
# # SIMULATION
# # -----------------------------------------------------------------------------
# # Configure simulation runs.
# run:
#
#   # Files output by such runs.
#   file:
#
#     # File caching the results of such runs.
#     cache: sim_init.betse
#
# # -----------------------------------------------------------------------------
# # PLOTS
# # -----------------------------------------------------------------------------
# # Configure simulation plotting.
# plot:
#
#   # Directories output by such plotting.
#   dir:
#
#     # Directory saving plotted images and movies.
#     media: sim_init.betse
#
# # -----------------------------------------------------------------------------
# # INTERNAL USE ONLY
# # -----------------------------------------------------------------------------
# # Avoid editing the following settings, which BETSE strictly requires for its
# # internal use only.
#
# # Configuration file version to which this file conforms. For reliable
# # comparability, this is stored as a string rather than float scalar.
# version: "0.0"
# '''
# '''
# Default YAML template for configuring simulations. Such string is formatted
# with  `{`- and `}`-delimited substrings intended to be expanded before use by a
# subsequent call to such string's format() method.
# '''

#FUXME: While useful, after some cotemplation we realize that we'd much rather
#have the default YAML content below saved to disk rather than embedded in a
#Python module. The former approach eliminates synchronization issues inherent
#to the latter approach (e.g., between changes Ally makes to local YAML files
#and their outdated contents in Python modules). So, that's clearly the road we
#should take. Unfortunately, doing so reliably raises several questions:
#
#* When installed with setuptools, how do we portably access such YAML file?
#* When installed with PyInstaller, how do we portably access such YAML file?
#
#Clearly, the codebase will need to transparently support both installation
#methods. To that end, let's get setuptools-based support working first and only
#worry about PyInstaller later.

# # Configuration file version to which this file conforms. For comparability,
# # this is split into major and minor components -- which *MUST* be integers.
# version:
#   major: 0
#   minor: 0

# -----------------------------------------------------------------------------
# EXTERNAL USE
# -----------------------------------------------------------------------------
# You are welcome to change any of the following settings.
# Version of the tissue simulation configuration format to which this file
# conforms.
# Tissue simulation configuration version to which this file adheres.
    # YAML Pathnames
    # ----------
    # For portability, such string embeds relative rather than absolute pathnames.
    # This instructs simulation objects elsewhere to save simulation output (e.g.,
    # results, plots) into such file's directory, permitting such directory to be
    # moved without breaking such configuration.
    # Get the default contents of YAML files configuring simulations.
# def get(
#     init_pickled_filename: str,
#     run_pickled_filename: str,
#     plot_output_dirname: str,
# ) -> str:
#     '''
#     Get a default YAML document for simulation configuration, with
#     the passed output configuration.
#     '''
#     assert isinstance(init_pickled_filename, str),\
#         '"{}" not a string.'.format(init_pickled_filename)
#     assert isinstance(run_pickled_filename, str),\
#         '"{}" not a string.'.format(run_pickled_filename)
#     assert isinstance(run_pickled_filename, str),\
#         '"{}" not a string.'.format(run_pickled_filename)
