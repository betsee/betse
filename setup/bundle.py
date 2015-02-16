#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
`betse`-specific `bundle` commands for `setuptools`.
'''

# ....................{ IMPORTS                            }....................
from setup.cmd import Command
# import os

# ....................{ COMMANDS                           }....................
def add_commands(setup_options: dict) -> None:
    '''
    Add `bundle` commands to the passed dictionary of `setuptools` options.
    '''
    assert isinstance(setup_options, dict),\
        '"{}" not a dictionary.'.format(setup_options)

    #FIXME: Common functionality. Contemplate a utility function. To implement
    #such function, we'd probably want such function to accept a list of class
    #objects rather than class names, and then simply access the "__name__"
    #attribute of each class object to dynamically obtain its name. Quite a bit
    #simpler than the current approach, when one considers it.

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

# ....................{ CLASSES                            }....................
#FIXME: Also make a "bundle_dir" class.

class bundle_file(Command):
    '''
    Create one platform-specific executable file in the top-level `dist`
    directory for each previously installed script.

    Each such file is created by running PyInstaller's external command
    `pyinstaller` with sane command-line arguments. Since PyInstaller does *not*
    currently (and probably never will) support cross-bundling, such files are
    formatted specific to and hence executable *only* under the currenty
    operating system. Specifically:

    * Under Linux, such files will be ELF (Executable and Linkable Format)
      binaries.
    * Under OS X, such files will be conventional ".app"-suffixed directories.
      (Of course, that's not a file. So sue us.)
    * Under Windows, such files will be conventional ".exe"-suffixed binaries.
    '''

    def finalize_options(self):
        '''
        Default undefined command-specific options to the options passed to the
        current parent command if any (e.g., `symlink`).
        '''
        # Copy the "install_dir" attribute from the existing "install_scripts"
        # attribute of a temporarily instantiated "symlink" object.
        #
        # Why? Because setuptools.
        self.set_undefined_options(
            'symlink', ('install_scripts', 'install_dir'))

    def run(self):
        '''Run the current command and all subcommands thereof.'''
        #FIXME: Run "pyinstaller."
        pass

# --------------------( WASTELANDS                         )--------------------
# from setuptools.command.install import install
# from setuptools.command.install_lib import install_lib
# from setuptools.command.install_scripts import install_scripts
# from distutils.errors import DistutilsFileError
    # Class Design
    # ----------
    # Despite subclassing the `install_scripts` class, this class does *not*
    # install scripts. This class subclasses such class merely to obtain access to
    # metadata on installed scripts (e.g., installation directory).

#FUXME: We may need to actually subclass "install" instead. No idea. Just try
#accessing "self.install_dir" below. If that fails, try "self.install_scripts".
#If that fails, try subclassing "install" instead and repeating such access
#attempts. Yes, this sucks. That's setuptools for you.
