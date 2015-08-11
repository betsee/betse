#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
`betse`-specific `symlink` commands for `setuptools`.
'''

# ....................{ IMPORTS                            }....................
from setup import util
from setuptools.command.install import install
from setuptools.command.install_lib import install_lib
from setuptools.command.install_scripts import install_scripts
from os import path
import os

# ....................{ COMMANDS                           }....................
def add_setup_commands(metadata: dict, setup_options: dict) -> None:
    '''
    Add `symlink` commands to the passed dictionary of `setuptools` options.
    '''
    util.add_setup_command_classes(
        metadata, setup_options,
        symlink, symlink_lib, symlink_scripts, unsymlink)

# ....................{ INSTALLERS                         }....................
class symlink(install):
    '''
    Editably install (e.g., in a symbolically linked manner) `betse` into the
    active Python 3 interpreter *without* performing dependency resolution.

    Unlike the default `develop` command, this command is suitable for
    system-wide installation.
    '''

    description =\
        'install a symlink rather than copy of this package (for development)'
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

    def run(self):
        '''Run the current command and all subcommands thereof.'''
        # If the current operating system is *NOT* POSIX-compatible, such system
        # does *NOT* provide conventional symbolic links. Raise an exception.
        util.die_if_os_non_posix()

        # Run all subcommands.
        for subcommand_name in self.get_sub_commands():
            self.run_command(subcommand_name)

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

    def run(self):
        util.die_if_os_non_posix()

        # Absolute path of betse's top-level Python package in the current
        # directory.
        package_dirname = path.join(os.getcwd(), self._setup_options['name'])

        # Absolute path of such symbolic link.
        symlink_filename = path.join(
            self.install_dir,
            self._setup_options['name'])

        # If such link currently exists, remove such link.
        if util.is_symlink(symlink_filename):
            util.remove_symlink(symlink_filename)

        # (Re)create such link.
        print('Symbolically linking "{}" to "{}".'.format(
            package_dirname, symlink_filename))
        os.symlink(package_dirname, symlink_filename)

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
        self.set_undefined_options('build', ('build_scripts', 'build_dir'))
        self.set_undefined_options(
            'symlink',
            ('install_scripts', 'install_dir'),
            ('force', 'force'),
            ('skip_build', 'skip_build'),
        )

# ....................{ UNINSTALLERS                       }....................
class unsymlink(install):
    '''
    Editably uninstall (e.g., in a symbolically linked manner) `betse` from
    the active Python 3 interpreter.
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
        self.install_lib_dir = None
        self.install_scripts_dir = None

    def finalize_options(self):
        '''
        Default undefined command-specific options to the options passed to the
        current parent command if any (e.g., `symlink`).
        '''
        # Copy attributes from a temporarily instantiated "symlink" object into
        # the current object under different attribute names.
        self.set_undefined_options(
            'symlink',
            ('install_lib',     'install_lib_dir'),
            ('install_scripts', 'install_scripts_dir'),
        )

    def run(self):
        '''Run the current command and all subcommands thereof.'''
        # If the current operating system is *NOT* POSIX-compatible, such system
        # does *NOT* provide conventional symbolic links. Raise an exception.
        util.die_if_os_non_posix()

        # Remove the installed library symbolic link.
        util.remove_symlink(path.join(
            self.install_lib_dir,
            self._setup_options['name'],
        ))

        # Remove all installed scripts.
        for script_basename, _, _ in util.command_entry_points(self):
            util.remove_file(path.join(
                self.install_scripts_dir, script_basename))

# --------------------( WASTELANDS                         )--------------------
        # Declare all attributes subsequently set by the call to
        # self.set_undefined_options(). (Such method raises exceptions if such
        # attributes have *NOT* been declared. This is rather terrible but par
        # for the setuptools course.)
        # Declare all attributes subsequently set by the call to
        # self.set_undefined_options(). (Such method raises exceptions if such
        # attributes have *NOT* been declared.)
    #FUXME: Insufficient. This obviously needs to uninstall *ALL* previously
    #symlinked scripts as well.

        # Absolute path of the library-specific symbolic link.
        # symlink_filename = path.join(
        #     self.install_lib_dir,
        #     self._setup_options['name'])
        #
        # # Remove such link.
        # remove_symlink(symlink_filename)
    # def finalize_options(self):
    #     '''
    #     Default undefined command-specific options to the options passed to the
    #     current parent command if any (e.g., `symlink`).
    #     '''
    #     # Copy the "install_dir" attribute from the existing "install_lib"
    #     # attribute of a temporarily instantiated "symlink" object.
    #     #
    #     # Welcome to setuptools hell.
    #     self.set_undefined_options(
    #         'symlink',
    #         ('install_lib', 'install_lib_dir'),
    #         ('install_scripts', 'install_scripts_dir'),
    #     )

    #FUXME: Insufficient. This obviously needs to uninstall *ALL* previously
    #symlinked scripts as well.
    # assert isinstance(setup_options, dict),\
    #     '"{}" not a dictionary.'.format(setup_options)
    #
    # # For the name of each command class to be registered as a new command...
    # for command_class_name in (
    #     'symlink', 'symlink_lib', 'symlink_scripts', 'unsymlink'):
    #     # Class object for the class with such name.
    #     command_class = globals()[command_class_name]
    #
    #     # Register such command.
    #     setup_options['cmdclass'][command_class_name] = command_class
    #
    #     # Expose the passed dictionary of "setuptools" options to such class by
    #     # adding a new private class field "_setup_options" to such class. While
    #     # merely passing such dictionary to instances of such classes would be
    #     # obviously preferable, setuptools and hence setuputils requires commands
    #     # be specified as uninstantiated classes rather than instances. Hence,
    #     # the current approach.
    #     command_class._setup_options = setup_options

    # user_options = [('install-dir=', 'd', 'directory to install to'),]
    # '''List of command-specific options.'''

    # def initialize_options(self):
    #     '''
    #     Initialize command-specific options.
    #     '''
    #     self.install_dir = None
# ....................{ GETTERS                            }....................
# def get_command_symlink_filename(setuptools_command):
#     '''
#     Get the absolute path of the symbolic link for `betse`'s top-level Python
#     package under the active Python 3 interpreter.
#
#     Parameters
#     ----------
#     setuptools_command : setuptools.Command
#         Command with which to inspect such interpreter.
#
#     Returns
#     ----------
#     string
#         Absolute path of such link for subsequent caller manipulation.
#     '''
#     return path.join(
#         setuptools_command.install_lib,
#         setuptools_command._setup_options['name'])

    # setup_options['cmdclass']['symlink'] = symlink
    # setup_options['cmdclass']['symlink_lib'] = symlink_lib
    # setup_options['cmdclass']['symlink_scripts'] = symlink_scripts
    # setup_options['cmdclass']['unsymlink'] = unsymlink

    # Expose such options to such classes by adding a new private class field
    # "_setup_options" to each such class. While merely passing such dictionary
    # to instances of such classes would be obviously preferable, setuptools and
    # hence setuputils requires commands be specified as uninstantiated classes
    # rather than instances. Hence, the current approach.
    # for command_class in setup_options['cmdclass']:
    #     command_class._setup_options = setup_options

# ....................{ GLOBALS                            }....................
# _setup_options = {}
# '''
# Dictionary of `setuptools` options, initialized by the call to `add_commands()`
# and subsequentled accessed from within command classes.
#
# While merely passing such dictionary to instances of such classes would be
# obviously preferable, `setuptools` requires uninstantiated classes rather than
# instances. Hence, the currently lackluster approach.
# '''
#     global _setup_options

    # Add the passed dictionary of `setuptools` options to .
    # for command_class in setup_options['cmdclass']:
    #FUXME: Rethink this. A class-based approach would probably be significantly
    #more intelligible, if slightly more verbose.

    # Define such functions as closures to provide such functions read-only
    # access to such dictionary.
