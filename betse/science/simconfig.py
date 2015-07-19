#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Simulation configuration in YAML format.
'''

#FIXME: Validate the versions of loaded configuration files.

# ....................{ IMPORTS                            }....................
from betse import pathtree
from betse.util.io import loggers
from betse.util.path import dirs, files, paths
import yaml

# ....................{ LOADERS                            }....................
def load(filename: str) -> dict:
    '''
    Load and validate the passed YAML file as a simulation configuration.

    After loading such file, this function validates the contents of such file.
    Specifically, this function raises an exception if:

    *
    '''
    # Contents of such file as a dictionary.
    config = None

    # Open the passed file for reading and read such file into a dictionary.
    with files.open_for_text_reading(filename) as yaml_file:
        config = yaml.load(yaml_file)

    #FIXME: Implement me.
    # Validate the contents of such file.

    # Get such dictionary.
    return config

# ....................{ WRITERS                            }....................
def write_default(filename: str) -> None:
    '''
    Write a default simulation configuration file in YAML format to the passed
    path _and_ copy all external resources (e.g., images) referenced by such
    configuration to the parent directory of such path.

    The resulting configuration will be usable as is with all high-level `betse`
    functionality requiring a valid configuration file (e.g., `betse world`).
    '''
    # Dirname, basename, and filetype of such file.
    dirname = paths.get_dirname(filename)
    basename = paths.get_basename(filename)
    filetype = paths.get_filetype(basename)

    # Log such creation.
    loggers.log_info(
        'Writing default simulation configuration to "{}".'.format(basename))

    # If such file already exists, fail. (For safety, we avoid silently
    # overwriting existing files.)
    files.die_if_file(filename)

    # If such filename is *NOT* suffixed by ".yaml", log a warning.
    if filetype != 'yaml':
        loggers.log_warning(
            'File "{}" filetype "{}" not "yaml".'.format(
                basename, filetype))

    # Create such file's parent directory, if needed.
    dirs.make_unless_dir(dirname)

    # Write the default configuration to such file.
    files.copy(pathtree.CONFIG_DEFAULT_FILENAME, filename)

    # Copy all external files referenced by such file to its parent directory.
    dirs.copy_into_target_dir(pathtree.DATA_GEOMETRY_DIRNAME, dirname)

# --------------------( WASTELANDS                         )--------------------
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

