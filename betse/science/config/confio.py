#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Simulation configuration file input and output (I/O) facilities.
'''

#FIXME: Validate the versions of loaded configuration files.

# ....................{ IMPORTS                            }....................
from betse import pathtree
from betse.util.io.log import logs
from betse.util.path import dirs, files, pathnames
from betse.util.path.dirs import DirOverwritePolicy
from betse.util.type.types import type_check  #, GeneratorType

# ....................{ WRITERS                            }....................
#FIXME: It should be feasible to replace this entire function (and hence remove
#this entire submodule) by:
#
#* Refactor the YamlFileABC.save() method signature to resemble this method's
#  signature as follows:
#
#    @type_check
#    def save(
#        # Mandatory parameters.
#        conf_filename: str,
#
#        # Optional parameters.
#        is_conf_file_overwritable: bool = False,
#        conf_subdir_overwrite_policy: DirOverwritePolicy = (
#            DirOverwritePolicy.SKIP_WITH_WARNING),
#    ) -> None:
#
#  Note our use of the "IGNORE_EXISTING" policy, which seems quite sensible,
#  safey, and sanitary for all YAML saving purposes.
#
#  *CAUTION*. Since the existing save() method performs overwriting by default,
#  we'll need to grep all calls to this method and explicitly pass the desired
#  parameters. Presumably, BETSEE already performs at least one such call.
#* Instantiating a "Parameters" object as follows:
#    p = Parameters().load(
#        conf_filename=pathtree.get_sim_config_default_filename())
#* Calling the save() method of this object.
#
#That's pretty obvious, frankly. Tragic that we didn't concoct it until now.
@type_check
def write_default(
    # Mandatory parameters.
    conf_filename: str,

    # Optional parameters.
    is_conf_overwritable: bool = False,
    is_data_overwritable: bool = False,
    # is_data_preservable: bool = False,
) -> None:
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
    is_conf_overwritable : optional[bool]
        ``True`` if an existing target YAML file is to be silently overwritten
        *or* ``False`` if an exception is to be raised if this file already
        exists. Defaults to ``False``.
    is_data_overwritable : optional[bool]
        ``True`` if existing target resources required by this target YAML file
        are to be silently overwritten *or* ``False`` if an exception is to be
        raised if any such resource already exists. Defaults to ``False``.
    is_data_preservable : optional[bool]
        ``True`` if existing target resources required by this target YAML file
        are to be preserved "as is" with a non-fatal warning *or* ``False`` if
        an exception is to be raised if any such resource already exists.
        Defaults to ``False``. The ``is_data_overwritable`` parameter takes
        precedence over this parameter. Specifically, if the
        ``is_data_overwritable`` parameter is ``True`` and this parameter is:
        * ``True``, an exception is raised. These two parameters conflict and
          hence *cannot* both be concurrently enabled.
        * ``False``, this parameter is silently ignored.

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
        is_overwritable=is_conf_overwritable,
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

            # Ignore (i.e., skip) each target resource of this subdirectory that
            # already exists with a non-fatal warning. This is purely a caller
            # convenience, permitting multiple configuration files with
            # different basenames to be created in the same parent directory
            # *WITHOUT* either raising exceptions on or silently overwriting
            # duplicate target resources.
            overwrite_policy=DirOverwritePolicy.SKIP_WITH_WARNING,

            # Ignore all empty ".gitignore" files in all subdirectories of this
            # source directory. These files serve as placeholders instructing
            # Git to track their parent subdirectories, but otherwise serve no
            # purpose. Preserving these files only invites end user confusion.
            ignore_basename_globs=('.gitignore',),
        )
