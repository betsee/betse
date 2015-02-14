#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Default YAML for simulation configuration.
'''

#FIXME: While useful, after some cotemplation we realize that we'd much rather
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
#FIXME: Validate the versions of loaded configuration files.

# ....................{ IMPORTS                            }....................
from betse import metadata

# ....................{ CONSTANTS ~ private                }....................
_TEMPLATE = '''
%YAML 1.1
---
# Default tissue simulation configuration.
#
# You are welcome to change any of the following settings.
#
# Pathnames
# ----------
# For portability, files and directories are configured below as relative rather
# than absolute pathnames (e.g., as "sim_init.betse" rather than
# "/Users/iamawesome/my_sim/sim_init.betse"). Relative pathnames instruct BETSE
# to save simulation results and plots relative to the directory containing this
# configuration file, permitting such directory to be moved elsewhere on the
# filesystem without "breaking" this configuration.
#
# YAML
# ----------
# As the "%YAML"-prefixed line above suggests, this file is formatted according
# to the YAML ("[Y]AML [A]in't [M]arkup [L]anguage") standard. YAML is a human-
# readable data serialization format appropriate for configuration files -- like
# this one. For further details on YAML syntax and semantics, see:
# https://en.wikipedia.org/wiki/YAML

# -----------------------------------------------------------------------------
# INITIALIZATION
# -----------------------------------------------------------------------------
# The following settings configure simulation initialization.
init:

  # Files output by simulation initialization.
  file:

    # Name of the file output by simulation initialization and subsequently
    # read by simulation runs.
    cache: sim_init.betse

# -----------------------------------------------------------------------------
# SIMULATION
# -----------------------------------------------------------------------------
# The following settings configure simulation runs.
run:

  # Files output by simulation runs.
  file:

    # Name of the file output by simulation runs and subsequently read by
    # simulation plotting.
    cache: sim_init.betse

# -----------------------------------------------------------------------------
# PLOTS
# -----------------------------------------------------------------------------
# The following settings configure simulation plotting.
plot:

  # Directories output by such plotting.
  dir:

    # Name of the directory to which such plotting outputs media files (e.g.,
    # images, animations).
    media: sim_init.betse

# -----------------------------------------------------------------------------
# INTERNAL USE ONLY
# -----------------------------------------------------------------------------
# Avoid editing the following settings, which BETSE strictly requires for its
# internal use only.

# Configuration file version to which this file conforms. For reliable
# comparability, this is stored as a string rather than float scalar.
version: "0.0"
'''
'''
Default YAML template for configuring simulations. Such string is formatted
with  `{`- and `}`-delimited substrings intended to be expanded before use by a
subsequent call to such string's format() method.
'''

# ....................{ GETTERS                            }....................
def get() -> str:
    '''
    Get the default YAML string for configuring simulations.

    Such string is intended to be subsequently serialized to disk (e.g., as a
    file suffixed by `.yaml`).
    '''
    return _TEMPLATE.format(
        program_name = metadata.NAME,
    )

# --------------------( WASTELANDS                         )--------------------
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

