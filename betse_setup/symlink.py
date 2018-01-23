#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
BETSE-specific `symlink` subcommands for `setuptools`.

## Microsoft Windows

Microsoft Windows does _not_ comply with POSIX standards and hence does _not_
generally support symbolic links.

While post-Vista versions of Microsoft Windows _do_ purport to support symbolic
links, the Windows version of the (Ana|Mini)conda Python distribution does _not_
appear to (at least, not reliably). Since this renders symbolic links useless
for standard Windows use, this module assumes Windows to _never_ support
symbolic links regardless of version.

Under Microsoft Windows, this module "fakes" symbolic link-based installation by
artifically prepending the Python-specific `sys.path` list of search dirnames
with the absolute path of the parent directory containing the top-level Python
package -- which largely has the same effect, albeit less resiliently. While
`sys.path` manipulation is (justifiably) frowned upon, no alternatives exist.
'''

# ....................{ IMPORTS                            }....................
import os
from betse_setup import buputil
from setuptools.command.install import install
from setuptools.command.install_lib import install_lib
from setuptools.command.install_scripts import install_scripts
from os import path

# ....................{ COMMANDS                           }....................
def add_setup_commands(metadata: dict, setup_options: dict) -> None:
    '''
    Add `symlink` subcommands to the passed dictionary of `setuptools` options.
    '''

    buputil.add_setup_command_classes(
        metadata, setup_options,
        symlink, symlink_lib, symlink_scripts, unsymlink)

# ....................{ CLASSES ~ install                  }....................
class symlink(install):
    '''
    Editably install (e.g., in a symbolically linked manner) `betse` into the
    active Python 3 interpreter *without* performing dependency resolution.

    Unlike the default `develop` command, this command is suitable for
    system-wide installation.
    '''

    # ..................{ ATTRIBUTES                         }..................
    description = (
        'install a symlink rather than copy of this package (for development)')
    '''
    Command description printed when running `./setup.py --help-commands`.
    '''


    sub_commands = [
        ('symlink_lib', None),
        ('symlink_scripts', None),
    ]
    '''
    Dictionary mapping command names to either:

    * A predicate returning a boolean indicating whether such command should be
      run under this run of the current command.
    * `None` signifying that such command should always be run.
    '''

    # ..................{ SUPERCLASS                         }..................
    def finalize_options(self):
        '''
        Default undefined command-specific options to the options passed to the
        current parent command if any (e.g., `install`).
        '''

        super().finalize_options()

        #FIXME: Replicate this functionality for the "install" command as well.

        # If the current system is OS X *AND* the OS X-specific Homebrew package
        # manager is installed...
        if buputil.is_os_os_x() and buputil.is_pathable('brew'):
             # Absolute path of Homebrew's top-level system-wide cellar
             # directory (e.g., "/usr/local/Cellar").
             brew_cellar_dir = buputil.get_command_output('brew', '--cellar')
             #print('Here!')

             # Absolute path of Homebrew's top-level system-wide directory
             # (e.g., "/usr/local").
             brew_dir = buputil.get_command_output('brew', '--prefix')

             # Absolute path of Homebrew's top-level system-wide binary
             # directory (e.g., "/usr/local/bin").
             brew_binary_dir = path.join(brew_dir, 'bin')

             # If this directory does not exist, raise an exception.
             buputil.die_unless_dir(brew_binary_dir)

             # If the directory to which wrappers will be installed is a Python-
             # specific subdirectory of this cellar directory (e.g.,
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
        '''Run the current command and all subcommands thereof.'''
        # If the current operating system is POSIX-incompatible, this system
        # does *NOT* support conventional symbolic links. See details above.
        if not buputil.is_os_posix():
            # Avoid circular import dependencies.
            from betse_setup import build

            # Print a non-fatal warning.
            buputil.output_warning(
                'Symbolic links require POSIX compatibility. '
                'Since the current platform is\n'
                'POSIX-incompatible (e.g., Windows), '
                'symbolic links will be faked with black magic.'
            )

            # Absolute path of the parent directory containing the top-level
            # "betse" package.
            parent_dirname = buputil.get_project_dirname()
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
            """.format(repr(parent_dirname)) + build.SCRIPT_TEMPLATE

        # Run all subcommands.
        for subcommand_name in self.get_sub_commands():
            self.run_command(subcommand_name)

# ....................{ CLASSES ~ install : subcommands    }....................
class symlink_lib(install_lib):
    '''
    Install the symbolic link for `betse`'s current editable installation,
    usually to a system-wide `site-packages` directory for the active Python 3
    interpreter.
    '''

    description = "install a symlink to this package's top-level module"
    '''
    Command description printed when running `./setup.py --help-commands`.
    '''


    def finalize_options(self):
        '''
        Default undefined command-specific options to the options passed to the
        current parent command if any (e.g., `symlink`).
        '''
        # Copy attributes from a temporarily instantiated "symlink" object into
        # the current object under different attribute names.
        self.set_undefined_options(
            'symlink', ('install_lib', 'install_dir'))

        # Default all remaining options.
        super().finalize_options()


    def run(self):
        # If the current operating system is POSIX-incompatible, such system
        # does *NOT* support conventional symbolic links. Return immediately.
        if not buputil.is_os_posix():
            return

        # Absolute path of betse's top-level Python package in the current
        # directory.
        package_dirname = path.join(os.getcwd(), self._setup_options['name'])

        # Absolute path of such symbolic link.
        symlink_filename = path.join(
            self.install_dir,
            self._setup_options['name'])

        # (Re)create such link.
        buputil.make_symlink(package_dirname, symlink_filename)


class symlink_scripts(install_scripts):
    '''
    Install all scripts wrapping `betse`'s current editable installation,
    usually to a system-wide directory in the current `${PATH}`.
    '''

    description =\
        'install scripts running this package without dependency checks'
    '''
    Command description printed when running `./setup.py --help-commands`.
    '''


    def finalize_options(self):
        '''
        Default undefined command-specific options to the options passed to the
        current parent command if any (e.g., `symlink`).
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

# ....................{ UNINSTALLERS                       }....................
class unsymlink(install):
    '''
    Editably uninstall (e.g., in a symbolically linked manner) `betse` from
    the active Python 3 interpreter.

    Attributes
    ----------
    install_package_dirname : str
        Absolute path of the directory to which our Python codebase was
        previously installed.
    install_wrapper_dirname : str
        Absolute path of the directory to which our wrapper scripts were
        previously installed.
    '''

    description =\
        'uninstall all installed symbolic links and scripts for this package'
    '''
    Command description printed when running `./setup.py --help-commands`.
    '''


    def initialize_options(self):
        '''
        Declare option-specific attributes subsequently initialized by
        `finalize_options()`.

        If this function is *not* defined, the default implementation of this
        method raises an inscrutable `distutils` exception. If such attributes
        are *not* declared, the subsequent call to
        `self.set_undefined_options()` raises an inscrutable `setuptools`
        exception. (This is terrible. So much hate.)
        '''
        super().initialize_options()
        self.install_package_dirname = None
        self.install_wrapper_dirname = None


    def finalize_options(self):
        '''
        Default undefined command-specific options to the options passed to the
        current parent command if any (e.g., `symlink`).
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
        # If the current operating system is POSIX-compatible, such system
        # supports symbolic links. In such case, remove the previously installed
        # symbolic link.
        if buputil.is_os_posix():
            buputil.remove_symlink(path.join(
                self.install_package_dirname,
                self._setup_options['name'],
            ))

        # Remove all installed scripts.
        for script_basename, _, _ in buputil.command_entry_points(self):
            buputil.remove_file(path.join(
                self.install_wrapper_dirname, script_basename))
