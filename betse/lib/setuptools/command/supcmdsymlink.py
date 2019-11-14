#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level custom ``symlink`` :mod:`setuptools` subcommands.

## Microsoft Windows

Microsoft Windows does *not* comply with POSIX standards and hence does *not*
generally support symbolic links.

While post-Vista versions of Microsoft Windows purport to support symbolic
links, the Windows version of the (Ana|Mini)conda Python distribution does
*not* appear to (at least, not reliably). Since this renders symbolic links
useless for standard Windows use, this module assumes Windows to *never*
support symbolic links regardless of operating system version.

Under Microsoft Windows, this module "fakes" symbolic link-based installation
by artifically prepending the Python-specific :attr:`sys.path` list of search
dirnames with the absolute path of the parent directory containing the
top-level Python package -- which largely has the same effect, albeit less
resiliently. While :attr:`sys.path` manipulation is (justifiably) frowned upon,
no plausible alternatives exist.
'''

# ....................{ IMPORTS                           }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To raise human-readable exceptions on missing mandatory
# dependencies, the top-level of this module may import *ONLY* from packages
# guaranteed to exist at installation time -- which typically means *ONLY*
# BETSE packages and stock Python packages.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

from betse.lib.setuptools.command import supcommand
from setuptools.command.install import install
from setuptools.command.install_lib import install_lib
from setuptools.command.install_scripts import install_scripts

# ....................{ ADDERS                            }....................
def add_subcommand(setup_options: dict, custom_metadata: dict) -> None:
    '''
    Add custom ``symlink`` :mod:`setuptools` subcommands to the passed
    dictionaries of :mod:`setuptools` options and arbirtrary metadata.
    '''

    # Who is Number One?
    supcommand.add_subcommand(
        setup_options, custom_metadata,
        symlink, symlink_lib, symlink_scripts, unsymlink)

# ....................{ SUBCOMMANDS                       }....................
class symlink(install):
    '''
    Editably install (e.g., in a symbolically linked manner) this application
    for the active Python interpreter *without* unnecessary (and occasionally
    harmful) dependency resolution.

    Unlike the default ``develop`` command, this command is suitable for
    system-wide installation.
    '''

    # ..................{ ATTRIBUTES                        }..................
    description = (
        'install a symlink rather than copy of this package (for development)')
    '''
    Command description printed when running ``./setup.py --help-commands``.
    '''


    sub_commands = [
        ('symlink_lib', None),
        ('symlink_scripts', None),
    ]
    '''
    Dictionary mapping command names to either:

    * A predicate returning a boolean indicating whether such command should be
      run under this run of the current subcommand.
    * ``None``, signifying that this subcommand should always be run.
    '''

    # ..................{ SUPERCLASS                        }..................
    def finalize_options(self):
        '''
        Default undefined command-specific options to the options passed to the
        current parent command if any (e.g., ``install``).
        '''

        # Defer heavyweight imports.
        from betse.util.os.brand import macos
        from betse.util.path import dirs, pathnames
        from betse.util.os.command import cmdrun, cmds

        # Finalize superclass options.
        super().finalize_options()

        #FIXME: Replicate this functionality for the "install" command as well.

        # If the current system is OS X *AND* the OS X-specific Homebrew package
        # manager is installed...
        if macos.is_macos() and cmds.is_cmd('brew'):
            # Absolute dirname of Homebrew's top-level system-wide cellar
            # directory (e.g., "/usr/local/Cellar").
            brew_cellar_dir = cmdrun.get_output_interleaved_or_die(
                command_words=('brew', '--cellar'))
            #print('Here!')

            # Absolute dirname of Homebrew's top-level system-wide directory
            # (e.g., "/usr/local").
            brew_dir = cmdrun.get_output_interleaved_or_die(
                command_words=('brew', '--prefix'))

            # Absolute dirname of Homebrew's top-level system-wide binary
            # directory (e.g., "/usr/local/bin").
            brew_binary_dir = pathnames.join(brew_dir, 'bin')

            # If this directory does not exist, raise an exception.
            dirs.die_unless_dir(brew_binary_dir)

            # If the directory to which wrappers will be installed is a
            # Python-specific subdirectory of this cellar directory (e.g.,
            # "/usr/local/Cellar/python3/3.5.0/Frameworks/Python.framework/Versions/3.5/bin"),
            # that subdirectory is unlikely to reside in the current ${PATH},
            # in which case wrappers installed to that subdirectory will remain
            # inaccessible. Correct this by forcing wrappers to be installed
            # to the Homebrew's conventional binary directory instead.
            if self.install_scripts.startswith(brew_cellar_dir):
                self.install_scripts = brew_binary_dir
                print('Detected Homebrew installation directory "{}".'.format(
                    brew_binary_dir))


    def run(self):
        '''
        Run the current command and all subcommands thereof.
        '''

        # Defer heavyweight imports.
        from betse.util.app.meta import appmetaone
        from betse.util.io import stderrs
        from betse.util.os.brand import posix

        # If the current operating system is POSIX-incompatible, this system
        # does *NOT* support conventional symbolic links. See details above.
        if not posix.is_posix():
            # Avoid circular import dependencies.
            from betse_setup import build

            # Print a non-fatal warning.
            stderrs.output_warning(
                'Symbolic links require POSIX compatibility. '
                'Since the current platform is\n'
                'POSIX-incompatible (e.g., Windows), '
                'symbolic links will be faked with black magic.'
            )

            # Absolute dirname of this application's project directory.
            project_dirname = appmetaone.get_app_meta().project_dirname
            # print('parent: ' + parent_dirname)

            # Prepend the template for subsequently installed entry points by a
            # Python statement "faking" symlink-based installation.
            build.SCRIPT_TEMPLATE = """
# The current operating system is POSIX-incompatible and hence does *NOT*
# support symlinks. To "fake" symlink-based installation, the standard list of
# search dirnames is prepended by the absolute path of the parent directory of
# the top-level "betse" package. For compatibility with third-party modules,
# the first entry of such list (i.e., the parent directory of this script) is
# preserved by inserting at index 1 rather than 0.
import sys
sys.path.insert(1, {})
            """.format(repr(project_dirname)) + build.SCRIPT_TEMPLATE

        # Run all subcommands.
        for subcommand_name in self.get_sub_commands():
            self.run_command(subcommand_name)

# ....................{ SUBCOMMANDS ~ subclasses          }....................
class symlink_lib(install_lib):
    '''
    Install the symbolic link for this application's current editable
    installation, usually to a system-wide ``site-packages`` directory for the
    active Python interpreter.
    '''

    description = "install a symlink to this package's top-level module"
    '''
    Command description printed when running ``./setup.py --help-commands``.
    '''


    def finalize_options(self):
        '''
        Default undefined command-specific options to the options passed to the
        current parent command if any (e.g., ``symlink``).
        '''

        # Copy attributes from a temporarily instantiated "symlink" object into
        # the current object under different attribute names.
        self.set_undefined_options(
            'symlink', ('install_lib', 'install_dir'))

        # Default all remaining options.
        super().finalize_options()


    def run(self):

        # Defer heavyweight imports.
        from betse.util.os.brand import posix
        from betse.util.os.shell import shelldir
        from betse.util.path import files, pathnames

        # If the current operating system is POSIX-incompatible, such system
        # does *NOT* support conventional symbolic links. Return immediately.
        if not posix.is_posix():
            return

        # Absolute path of betse's top-level Python package in the current
        # directory.
        package_dirname = pathnames.join(
            shelldir.get_cwd_dirname(), self._setup_options['name'])

        # Absolute path of such symbolic link.
        symlink_filename = pathnames.join(
            self.install_dir,
            self._setup_options['name'])

        #FIXME: Define a new files.make_symlink() function resembling the
        #existing buputils.make_symlink() function.
        # (Re)create such link.
        files.make_symlink(package_dirname, symlink_filename)


class symlink_scripts(install_scripts):
    '''
    Install all scripts wrapping this application's current editable
    installation, usually to a system-wide directory in the current
    ``${PATH}``.
    '''

    description = (
        'install scripts running this package without dependency checks')
    '''
    Command description printed when running ``./setup.py --help-commands``.
    '''


    def finalize_options(self):
        '''
        Default undefined command-specific options to the options passed to the
        current parent command if any (e.g., ``symlink``).
        '''

        # Copy attributes from a temporarily instantiated "symlink" object into
        # the current object under different attribute names.
        self.set_undefined_options(
            'build',
            ('build_scripts', 'build_dir'),
        )
        self.set_undefined_options(
            'symlink',
            ('install_scripts', 'install_dir'),
            ('force', 'force'),
            ('skip_build', 'skip_build'),
        )

        # Default all remaining options.
        super().finalize_options()

# ....................{ UNINSTALLERS                      }....................
class unsymlink(install):
    '''
    Editably uninstall (e.g., in a symbolically linked manner) this application
    from the active Python interpreter.

    Attributes
    ----------
    install_package_dirname : str
        Absolute path of the directory to which our Python codebase was
        previously installed.
    install_wrapper_dirname : str
        Absolute path of the directory to which our wrapper scripts were
        previously installed.
    '''

    description = (
        'uninstall all installed symbolic links and scripts for this package')
    '''
    Command description printed when running `./setup.py --help-commands`.
    '''


    def initialize_options(self):
        '''
        Declare option-specific attributes subsequently initialized by
        :meth:`finalize_options`.

        If this function is *not* defined, the default implementation of this
        method raises an inscrutable :mod:`distutils` exception. If these
        attributes are *not* declared, the subsequent call to
        :meth:`set_undefined_options` raises an inscrutable :mod:`setuptools`
        exception. (This is terrible. So much hate.)
        '''
        super().initialize_options()
        self.install_package_dirname = None
        self.install_wrapper_dirname = None


    def finalize_options(self):
        '''
        Default undefined command-specific options to the options passed to the
        current parent command if any (e.g., ``symlink``).
        '''

        # Copy attributes from a temporarily instantiated "symlink" object into
        # the current object under different attribute names.
        self.set_undefined_options(
            'symlink',
            ('install_lib',     'install_package_dirname'),
            ('install_scripts', 'install_wrapper_dirname'),
        )

        # Default all remaining options.
        super().finalize_options()


    def run(self):
        '''Run the current command and all subcommands thereof.'''

        # Defer heavyweight imports.
        from betse.util.os.brand import posix
        from betse.util.path import files, pathnames

        # If the current operating system is POSIX-compatible, such system
        # supports symbolic links. In such case, remove the previously
        # installed symbolic link.
        if posix.is_posix():
            #FIXME: Define a new files.remove_symlink() function resembling the
            #existing buputils.remove_symlink() function.
            files.remove_symlink(pathnames.join(
                self.install_package_dirname,
                self._setup_options['name'],
            ))

        # Remove all installed scripts.
        for script_basename, _, _ in (
            supcommand.iter_command_entry_points(self)):
            files.remove_file(pathnames.join(
                self.install_wrapper_dirname, script_basename))
