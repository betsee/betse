#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Utility functions converting the default POSIX-specific YAML simulation
configuration bundled with BETSE into a platform-specific YAML simulation
configuration appropriate for use by the current user on the current system.
'''

# ....................{ IMPORTS                            }....................
from betse import pathtree

#FIXME: This isn't particularly good. Contemplate alternatives, please.
#Actually, this will probably go away when we refactor to use SimConfigWrapper.
from betse.science.config.confio import _write_check

from betse.util.io.log import logs
from betse.util.path import paths, dirs, files
from betse.util.type.types import type_check

# ....................{ WRITERS                            }....................
#FIXME: Refactor this function as follows (in order):
#
#* Refactor this function to internally create and return a "SimConfigWrapper"
#  instance. Since this function is only called once elsewhere in the codebase,
#  this should be relatively simple. Do so now, however, before matters become
#  more... entangled.
@type_check
def write(config_filename: str) -> None:
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
    BetseFileException
        If this file already exists.
    '''

    # Announce the ugly shape of things to come.
    logs.log_info('Writing default simulation configuration...')

    # Validate this file *BEFORE* creating this file's parent directory if
    # needed *BEFORE* creating this file.
    _write_check(config_filename)
    _write_dir(config_filename)
    _write_file(config_filename)


@type_check
def _write_dir(config_filename: str) -> None:
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

    # Parent directory of this file if any or the current directory otherwise.
    target_dirname = paths.get_dirname_or_current_dirname(config_filename)

    # Create this directory if needed.
    dirs.make_unless_dir(target_dirname)

    # For the absolute path of each source subdirectory containing assets
    # required by this file...
    for source_asset_dirname in pathtree.get_data_default_asset_dirnames():
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


@type_check
def _write_file(config_filename: str) -> None:
    '''
    Write a default YAML simulation configuration file to the passed path,
    assumed to _not_ already exist.

    Parameters
    ----------
    config_filename : str
        Absolute or relative path of the YAML file to be written.
    '''

    #FIXME: Ideally, we should be using ruamel.yaml to munge YAML data in a
    #well-structured and hence sane manner rather than the admittedly crude
    #"sed"-like approach leveraged below. Unfortunately, given the complex
    #nature of this data, it's unclear whether or not ruamel.yaml would
    #adequately preserve the entirety of this data in a roundtrip manner. For
    #now, the "sed"-like approach prevails.

    # Write the default configuration to this file, modifying the latter with
    # "sed"-like global string substitution as detailed above.
    files.replace_substrs(
        filename_source=pathtree.get_sim_config_default_filename(),
        filename_target=config_filename,
        replacements=(
            # Prevent static plots from being displayed by default.
            (r'^(\s*plot after solving:\s+)True\b(.*)$', r'\1False\2'),
        ),
    )
