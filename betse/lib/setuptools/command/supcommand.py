#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level :mod:`setuptools` command facilities.
'''

# ....................{ IMPORTS                           }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To raise human-readable exceptions on missing mandatory
# dependencies, the top-level of this module may import *ONLY* from packages
# guaranteed to exist at installation time -- which typically means *ONLY*
# BETSE packages and stock Python packages.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

from betse.util.path import pathnames
from betse.util.type.types import type_check, GeneratorType, MappingType
from pkg_resources import Distribution, PathMetadata
from setuptools import Command
from setuptools.command.develop import VersionlessRequirement

# ....................{ TYPES                             }....................
SetuptoolsCommandDistributionTypes = (Distribution, VersionlessRequirement)
'''
**Distribution** (i.e., high-level object encapsulating metadata for a
:mod:`setuptools`-installed Python project) as commonly passed to methods
called by the :mod:`setuptools.command.easy_install.easy_install` class.

Specifically, if the end user invoked the :mod:`setuptools` subcommand:

* ``develop``, this object is an instance of the
  :mod:`setuptools`-specific :class:`VersionlessRequirement` class.
  Confusingly, note that this class wraps the underlying
  :mod:`pkg_resources`-specific :class:`Distribution` class as a class proxy
  transparently stripping versioning from this distribution's name (e.g., by
  truncating ``foo==1.0`` to merely ``foo``).
* ``install``, this object is an instance of the :mod:`pkg_resources`-specific
  :class:`Distribution` class. Confusingly, note that this class has no
  relation whatsoever to the identically named
  :class:`distutils.dist.Distribution` and
  :class:`setuptools.dist.Distribution` classes.

Why, :mod:`setuptools:`. Why.
'''

# ....................{ ADDERS                            }....................
@type_check
def add_subcommand(
    setup_options: MappingType,
    custom_metadata: MappingType,

    #FIXME: Type-check each such command to be a "Command". For unknown
    #reasons, Python is currently complaining that "freeze_dir" is not a
    #"Command", despite clearly being a "Command". </collective_shrug>
    *subcommands
) -> None:
    '''
    Define one custom :mod:`setuptools` subcommand for each passed class,
    configured by the passed dictionaries of :mod:`setuptools` options and
    arbirtrary metadata.

    For simplicity, the name of each such subcommand will be the name of the
    corresponding class. Hence, the names of these classes are recommended to
    be terse lowercase strings (e.g., ``freeze``, ``symlink``).

    Parameters
    ----------
    setup_options : dict
        Dictionary of **:mod:`setuptools` options** (i.e., :mod:`setuptools`-
        rather than application-specific), mapping from the names of parameters
        accepted by the standard :func:`setuptools.setup` function to the
        values of those parameters. For each passed command subclass, this
        function adds a new entry to the ``cmdclass`` key of this dictionary
        exposing that subclass to :mod:`setuptools`.
    custom_metadata : dict
        Dictionary of **arbitrary metadata** (i.e., application- rather than
        :mod:`setuptools`-specific), mapping from arbitrary keys to values.
        This metadata is intended to inform custom subcommands (e.g.,
        ``freeze_file``) of supplementary metadata *not* already declared by
        the :data:`setup_options` dictionary. Since :mod:`setuptools` raises
        fatal exceptions on detecting unrecognized keys in the passed
        ``setup_options`` dictionary, unrecognized keys are added to this
        dictionary instead.
    subcommands : tuple
        Tuple of zero or more subclasses of the standard
        :class:`setuptools.Command` superclass to be defined, each of which is
        assumed to implement an application-specific :mod:`setuptools`
        subcommand runnable by end users from the command line (e.g., by
        passing any subcommand name to the top-level ``setup.py`` script).
    '''

    # If this is the first set of custom commands to be defined and hence the
    # first call to this function, initialize the corresponding
    # setuptools.setup() option to the empty dictionary.
    if 'cmdclass' not in setup_options:
        setup_options['cmdclass'] = {}

    # For each such command class...
    for subcommand in subcommands:
        # Add this command class as a new command of the same name.
        setup_options['cmdclass'][subcommand.__name__] = subcommand

        # Expose the passed dictionaries to this class by monkey-patching
        # application-specific private class variables into these classes.
        # While passing these dictionaries to instances of this class (e.g., on
        # instantiation) would be ideal, distutils and hence setuptools
        # requires commands to be registered as classes rather than instances.
        subcommand._custom_metadata = custom_metadata
        subcommand._setup_options = setup_options

# ....................{ ITERATORS                         }....................
@type_check
def iter_subcommand_entry_points(subcommand: Command) -> GeneratorType:
    '''
    Generator yielding a 3-tuple detailing each wrapper script installed for
    the distribution described by the passed :mod:`setuptools` subcommand.

    See Also
    ----------
    :func:`iter_package_distribution_entry_points`
        Further details.
    '''

    # Make a "pkg_resources"-specific distribution from the passed command.
    # Yes, this code was ripped wholesale from the run() method defined by
    # module "setuptools.command.install_scripts". Yes, we don't know how it
    # works. "Frankly, Mam, we don't give a damn."
    #
    # It should be noted that all commands have an attribute "distribution".
    # Naturally, this is a setuptools-specific distribution that has literally
    # *NOTHING* to do with pkg_resources-style distributions.
    #
    # Die, setuptools. Die!
    ei_cmd = subcommand.get_finalized_command('egg_info')
    distribution = Distribution(
        ei_cmd.egg_base,
        PathMetadata(ei_cmd.egg_base, ei_cmd.egg_info),
        ei_cmd.egg_name,
        ei_cmd.egg_version,
    )

    # Defer to the generator provided by such function.
    yield from iter_package_distribution_entry_points(distribution)


@type_check
def iter_package_distribution_entry_points(
    distribution: SetuptoolsCommandDistributionTypes) -> GeneratorType:
    '''
    Generator iteratively yielding a 3-tuple describing each wrapper script
    installed for the passed distribution.

    Parameters
    ----------
    distribution : SetuptoolsCommandDistributionTypes
        **Distribution** (i.e., high-level object encapsulating metadata for a
        :mod:`setuptools`-installed Python project).

    Yields
    ----------
    (str, str, EntryPoint)
        3-tuple ``(script_basename, ui_type, entry_point)`` such that:

        * ``script_basename` is this script's basename (e.g., ``betse``). To
          simplify integration with the downstream setuptools API (e.g., the
          :meth:`setuptools.command.easy_install.ScriptWriter.get_script_args`
          method), this basename is typically *not* suffixed by a
          platform-specific filetype (e.g., ``.exe`` under vanilla or Cygwin
          Microsoft Windows).
        * ``ui_type`` is this script's interface type string, guaranteed to be
          either:

          * If this script is console-specific, ``console``.
          * Else, ``gui``.

        * ``entry_point`` is this script's :class:`pkg_resources.EntryPoint`
          object, whose attributes specify the module to be imported and
          function to be run by this script.
    '''

    # For each type of script wrapper...
    for script_type in 'console', 'gui':
        script_type_group = script_type + '_scripts'

        # For each script of this type...
        for script_basename, entry_point in (
            distribution.get_entry_map(script_type_group).items()):
            # If this basename is *NOT* a basename, raise an exception. Note
            # that similar validation is performed by the
            # ScriptWriter.get_args() class method inspiring this function.
            pathnames.die_unless_basename(script_basename)

            # Yield this 3-tuple. To simplify integration with the downstream
            # setuptools API, do *NOT* sanitize_snakecase this script's
            # basename by calling sanitize_command_basename(). Since that API
            # already implicitly suffixes this basename by ".exe", doing so
            # here would erroneously result in this basename being suffixed by
            # ".exe.exe".
            yield script_basename, script_type, entry_point
