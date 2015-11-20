#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

#FIXME: Executables output under OS X and Windows pretty much *MUST* be signed.
#This looks to be fairly trivial under Windows. OS X, however, is another kettle
#of hideous fish. In any case, everyone else has already solved this, so we just
#need to leverage the following detailed recipes:
#
#* https://github.com/pyinstaller/pyinstaller/wiki/Recipe-Win-Code-Signing
#* https://github.com/pyinstaller/pyinstaller/wiki/Recipe-OSX-Code-Signing

#FIXME: Embed Windows-specific version metadata in such executables. This is a
#fairly bizarre process, which we've documented in "pyinstaller.yaml". It's
#hardly crucial for now, but will be important at some point.
#FIXME: Embed OS X-specific version metadata in such executables. We'll want to
#detect whether the current OS is OS X and, if so, manually overwrite the
#autogenerated "myapp.app/Contents/Info.plist" file with one of our own
#devising. Not terribly arduous... in theory.

#FIXME: Embed ".ico"-suffixed icon files in such executables. PyInstaller
#provides simple CLI options for this; we simply need to create such icons.
#Contemplating mascots, how about the ever-contemplative BETSE cow?
#
#Sadly, the icon formats required by OS X and Windows appear to conflict.
#Windows icon files have filetype ".ico" (and appear to support only one
#embedded icon), whereas OS X icon files have filetype ".icns" (and appear to
#support multiple embedded icons). To compound matters, "pyinstaller" provides
#only one option "--icon" for both, probably implying that we'll need to
#dynamically detect whether the current system is OS X or Windows and respond
#accordingly (i.e., by passing the appropriate system-specific icon file).

#FIXME: Contribute back to the community. Contemplate a stackoverflow answer.
#(We believe we may have espied an unanswered question asking about query words
#"pyinstaller setuptools integration". Huzzah!) We should note that PyInstaller
#will probably be unable to find the imports of setuptools-installed scripts,
#due to the obfuscatory nature of such scripts. See the following for a
#reasonable solution:
#    https://github.com/pyinstaller/pyinstaller/wiki/Recipe-Setuptools-Entry-Point

'''
`betse`-specific `freeze` subcommands for `setuptools`.
'''

# ....................{ IMPORTS                            }....................
from abc import ABCMeta, abstractmethod
from distutils.errors import DistutilsExecError
from os import path
from pkg_resources import EntryPoint
from setup import util
from setuptools import Command
import os

# ....................{ COMMANDS                           }....................
def add_setup_commands(metadata: dict, setup_options: dict) -> None:
    '''
    Add `freeze` subcommands to the passed dictionary of `setuptools` options.
    '''
    util.add_setup_command_classes(
        metadata, setup_options, freeze_dir, freeze_file)

# ....................{ CLASSES                            }....................
class freeze(Command, metaclass = ABCMeta):
    '''
    Abstract command class creating either one platform-specific executable file
    *or* one directory containing such file in the top-level `dist` directory
    for each previously installed script of the current application.

    Each such file is created by running the external `pyinstaller` command with
    sane command-line arguments. Since such command does *not* support
    **cross-bundling** (i.e., creation of executable files for operating systems
    other than the current), such files apply *only* to the current such system.
    Specifically:

    * Under Linux, such files will be ELF (Executable and Linkable Format)
      binaries.
    * Under OS X, such files will be conventional ".app"-suffixed directories.
      (Of course, that's not a file. So sue us.)
    * Under Windows, such files will be conventional ".exe"-suffixed binaries.

    Attributes
    ----------
    clean : bool
        True if the user passed the `--clean` option to the current `setuptools`
        command.
    install_dir : str
        Absolute path of the directory to which our wrapper scripts were
        previously installed.
    _pyinstaller_command : list
        List of all shell words of the PyInstaller command to be run.
    _pyinstaller_spec_filename : str
        Relative path of the PyInstaller spec file converting the current
        platform-independent wrapper script into a platform-specific executable.
    _pyinstaller_dist_dirname : str
        Relative path of the PyInstaller directory to which the final executable
        (as either a file or directory) is written.
    _pyinstaller_hooks_dirname : str
        Relative path of the input hooks subdirectory.
    '''

    user_options = [
        ('clean', None,
         'clean PyInstaller cache of temporary paths before building'),
        ('debug', None,
         'print debug messages during PyInstaller bootloader startup'),
    ]
    '''
    List of 3-tuples specifying command-line options accepted by this command.

    For each such option, an attribute of the same name as such option's long
    form _must_ be explicitly initialized in the `initialize_options()` method.
    `setuptools` fails to recognize options for which this is _not_ the case.
    (You fail a simple sanity check yet again, `setuptools`.)

    See Also
    ----------
    http://ilostmynotes.blogspot.ca/2009/04/python-distutils-installer-and.html
        Inarguably, the best (albeit unofficial) documentation on such list.
    '''

    # ..................{ ATTRIBUTES                         }..................
    # Advice to be appended to exceptions raised below.
    EXCEPTION_ADVICE = (
        'Consider running either:\n'
        '\tsudo python3 setup.py install\n'
        '\tsudo python3 setup.py symlink')

    # ..................{ SUPERCLASS                         }..................
    def initialize_options(self):
        '''
        Declare option-specific attributes subsequently initialized by
        `finalize_options()`.

        If this function is _not_ defined, the default implementation of this
        method raises an inscrutable `distutils` exception. If such attributes
        are _not_ declared, the subsequent call to
        `self.set_undefined_options()` raises an inscrutable `setuptools`
        exception. (This is terrible. So much hate.)
        '''
        #FIXME: To circumvent PyInstaller caching issues, we currently force
        #"self.clean = True". Once PyInstaller caching is sufficiently reliable
        #to reasonably permit reuse of cached metadata, revert this back to
        #"self.clean = False".

        # Option-specific public attributes. For each option declared by the
        # "user_options" list above, a public attribute of the same name as such
        # option's long form *MUST* be initialized here to its default value.
        self.clean = True
        self.debug = False

        # setuptools-specific public attributes.
        self.install_dir = None

        # Custom private attributes.
        self._pyinstaller_command = None
        self._pyinstaller_spec_filename = None
        self._pyinstaller_dist_dirname = None
        self._pyinstaller_hooks_dirname = None

    def finalize_options(self):
        '''
        Default undefined command-specific options to the options passed to the
        current parent command if any (e.g., `symlink`).
        '''
        # Copy attributes from a temporarily instantiated "symlink" object into
        # the current object under different attribute names.
        self.set_undefined_options(
            'symlink', ('install_scripts', 'install_dir'))

    def run(self):
        '''Run the current command and all subcommands thereof.'''
        # List of all shell words of the PyInstaller command to be run.
        self._init_pyinstaller_command()

        # True if the current distribution has at least one entry point.
        is_entry_point = False

        # Freeze each previously installed script wrapper.
        for script_basename, script_type, entry_point in\
            util.command_entry_points(self):
            # Note at least one entry point to be installed.
            is_entry_point = True

            # Validate this wrapper's entry point.
            # freeze._check_entry_point(entry_point)

            # Relative path of the output frozen executable file or directory,
            # created by stripping the suffixing ".exe" filetype on Windows.
            frozen_pathname = path.join(
                self._pyinstaller_dist_dirname,
                util.get_path_sans_filetype(script_basename))

            # If cleaning and such path exists, remove such path *BEFORE*
	    # validating such path. Why? Such path could be an existing file
	    # and the current command freezing to a directory (or vice versa),
	    # in which case subsequent validation would raise an exception.
            if self.clean and util.is_path(frozen_pathname):
                util.remove_path(frozen_pathname)

            # Validate such path.
            self._check_frozen_path(frozen_pathname)

            # Set all environment variables used to communicate with the BETSE-
            # specific PyInstaller specification file run below.
            self._set_environment_variables(
                script_basename, script_type, entry_point)

            # Run the desired PyInstaller command.
            self._run_pyinstaller(
                script_basename, script_type, entry_point)

            # Report such results to the user.
            frozen_pathtype =\
                'directory' if util.is_dir(frozen_pathname) else 'file'
            print('Froze {} "{}".\n'.format(frozen_pathtype, frozen_pathname))

            #FIXME: Excise when beginning GUI work.
            break

        # If no entry points are registered for the current distribution, raise
        # an exception.
        if not is_entry_point:
            raise DistutilsExecError(
                'No entry points found. {}'.format(freeze.EXCEPTION_ADVICE))

    # ..................{ INITIALIZERS                       }..................
    def _init_pyinstaller_command(self) -> None:
        '''
        Initialize the list of all shell words of the PyInstaller command to be
        run.
        '''
        # Basename of the PyInstaller command to be run. To avoid confusion with
        # non-Windows executables in the current ${PATH} when running under Wine
        # emulation, accept only Windows executables when running under Windows.
        pyinstaller_command_name = 'pyinstaller'
        if util.is_os_windows():
            pyinstaller_command_name += '.exe'

        # If PyInstaller is not found, fail.
        util.die_unless_pathable(
            pyinstaller_command_name,
            'PyInstaller not installed or "{}" not in the current ${{PATH}}.'.format(
                pyinstaller_command_name))

        # List of all shell words of the PyInstaller command to be run.
        self._pyinstaller_command = [pyinstaller_command_name]

        # Relative path of the top-level PyInstaller directory.
        pyinstaller_dirname = 'freeze'

        # Relative path of the PyInstaller spec file converting such
        # platform-independent script into a platform-specific executable.
        self._pyinstaller_spec_filename = path.join(
            pyinstaller_dirname, '.spec')

        # Relative path of the final output subdirectory.
        self._pyinstaller_dist_dirname = path.join(pyinstaller_dirname, 'dist')

        # Relative path of the input hooks subdirectory.
        self._pyinstaller_hooks_dirname = path.join(pyinstaller_dirname, 'hooks')

        # Relative path of the intermediate build subdirectory.
        pyinstaller_work_dirname = path.join(pyinstaller_dirname, 'build')

        # Create such hooks subdirectory if not found, as failing to do so
        # will induce fatal PyInstaller errors.
        util.make_dir_unless_found(self._pyinstaller_hooks_dirname)

        # Append all PyInstaller command options common to running such command
        # for both reuse and regeneration of spec files. (Most such options are
        # specific to the latter only and hence omitted.)
        self._pyinstaller_command.extend([
            # Overwrite existing output paths under the "dist/" subdirectory
            # without confirmation, the default behaviour.
            '--noconfirm',

            # Non-default PyInstaller directories.
            '--workpath=' + util.shell_quote(pyinstaller_work_dirname),
            '--distpath=' + util.shell_quote(self._pyinstaller_dist_dirname),

            # Non-default log level.
            '--log-level=DEBUG',
        ])

        # Forward all custom boolean options passed by the user to the current
        # setuptools command (e.g., "--clean") to the "pyinstaller" command.
        if self.clean:
            self._pyinstaller_command.append('--clean')
        if self.debug:
            self._pyinstaller_command.extend((
                '--debug',

                # UPX-based compression uselessly consumes non-trivial time
                # (especially under Windows, where process creation is fairly
                # heavyweight) when freezing debug binaries. To optimize and
                # simplify debugging, such compression is disabled.
                '--noupx',
            ))
            util.output_warning(
                'Enabling bootloader debug messages.')
            util.output_warning(
                'Disabling UPX-based compression.')
        # If *NOT* debugging and UPX is *NOT* found, print a non-fatal warning.
        # While optional, freezing in the absence of UPX produces uncompressed
        # and hence considerably larger executables.
        elif not util.is_pathable('upx'):
            util.output_warning(
                'UPX not installed or "upx" not in the current ${PATH}.')
            util.output_warning(
                'Frozen binaries will *NOT* be compressed.')

    # ..................{ SETTERS                            }..................
    def _set_environment_variables(
        self, script_basename: str, script_type: str, entry_point: str) -> None:
        '''
        Set all environment variables used to communicate with the BETSE-
        specific PyInstaller specification file run in a separate process, given
        the passed arguments yielded by the `command_entry_points()` generator.

        While hardly ideal, PyInstaller appears to provide no other means of
        communicating with such file.
        '''
        # Absolute path of the entry module.
        #
        # This module's relative path to the top-level project directory is
        # obtained by converting the entry point specifier defined by "setup.py"
        # for the current entry point (e.g., "betse.gui.guicli:main") into a
        # platform-specific path. Sadly, setuptools provides no cross-platform
        # API for reliably obtaining the absolute path of the corresponding
        # script wrapper. Even if it did, such path would be of little use under
        # POSIX-incompatible platforms (e.g., Windows), where these wrappers are
        # binary blobs rather than valid Python scripts.
        #
        # Instead, we reverse-engineer the desired path via brute-force path
        # manipulation. Thus burns out another tawdry piece of my soul.
        module_filename = path.join(
            util.get_project_dirname(),
            entry_point.module_name.replace('.', path.sep) + '.py')

        # Ensure such module exists.
        util.die_unless_file(module_filename)

        # Such path.
        os.environ['__FREEZE_MODULE_FILENAME'] = module_filename

        # Whether to freeze in "one-file" or "one-directory" mode.
        os.environ['__FREEZE_MODE'] = self._get_freeze_mode()

        # Whether to freeze a CLI- or GUI-based application.
        os.environ['__FREEZE_INTERFACE_TYPE'] = script_type

    # ..................{ RUNNERS                            }..................
    def _run_pyinstaller(
        self, script_basename: str, script_type: str, entry_point: str) -> None:
        '''
        Run the desired PyInstaller command.
        '''
        # If such spec exists, instruct PyInstaller to reuse rather than
        # recreate such file, thus preserving edits to such file.
        if util.is_file(self._pyinstaller_spec_filename):
            print('Reusing spec file "{}".'.format(
                self._pyinstaller_spec_filename))

            # Append the relative path of such spec file.
            self._pyinstaller_command.append(
                util.shell_quote(self._pyinstaller_spec_filename))

            # Freeze such script with such spec file.
            util.die_unless_command_succeeds(*self._pyinstaller_command)
        # Else, instruct PyInstaller to (re)create such spec file.
        else:
            # Absolute path of the directory containing such files.
            pyinstaller_spec_dirname = util.get_path_dirname(
                self._pyinstaller_spec_filename)

            # Absolute path of the current script wrapper.
            script_filename = path.join(
                self.install_dir, script_basename)
            util.die_unless_file(
                script_filename, 'File "{}" not found. {}'.format(
                    script_filename, freeze.EXCEPTION_ADVICE))

            # Inform the user of such action *AFTER* the above validation.
            # Since specification files should typically be reused rather
            # than regenerated, do so as a non-fatal warning.
            util.output_warning(
                'Generating spec file "{}".'.format(
                    self._pyinstaller_spec_filename))

            # Append all options specific to spec file generation.
            self._pyinstaller_command.extend([
                # If this is a console script, configure standard input and
                # output for console handling; else, do *NOT* and, if the
                # current operating system is OS X, generate an ".app"-suffixed
                # application bundle rather than a customary executable.
                '--console' if script_type == 'console' else '--windowed',

                # Non-default PyInstaller directories.
                '--additional-hooks-dir=' + util.shell_quote(
                    self._pyinstaller_hooks_dirname),
                '--specpath=' + util.shell_quote(pyinstaller_spec_dirname),
            ])

            # Append all subclass-specific options.
            self._pyinstaller_command.extend(self._get_pyinstaller_options())

            # Append the absolute path of such script.
            self._pyinstaller_command.append(util.shell_quote(script_filename))

            # Freeze such script and generate a spec file.
            util.die_unless_command_succeeds(*self._pyinstaller_command)

            # Absolute path of such file.
            script_spec_filename = path.join(
                pyinstaller_spec_dirname, script_basename + '.spec')

            # Rename such file to have the basename expected by the prior
            # conditional on the next invocation of this setuptools command.
            #
            # Note that "pyinstaller" accepts an option "--name" permitting
            # the basename of such file to be specified prior to generating
            # such file. Unfortunately, such option *ALSO* specifies the
            # basename of the generated executable. While the former is reliably
            # renamable, the former is *NOT* (e.g., due to code signing). Hence,
            # such file is manually renamed without passing such option.
            util.move_file(
                script_spec_filename, self._pyinstaller_spec_filename)

    # ..................{ CLASS                              }..................
    @classmethod
    def _check_entry_point(entry_point: EntryPoint):
        '''
        Validate the passed entry point, describing the current script wrapper
        to be frozen.
        '''
        assert isinstance(entry_point, EntryPoint),\
            '"{}" not an entry point.'.format(entry_point)

        # If this entry module is unimportable, raise an exception.
        if not util.is_module(entry_point.module_name):
            raise ImportError(
                'Entry module "{}" unimportable. {}'.format(
                entry_point.module_name, freeze.EXCEPTION_ADVICE))

        # If this entry module's basename is *NOT* "__main__", print a
        # non-fatal warning. For unknown reasons, attempting to freeze
        # customary modules usually causes the frozen executable to reduce
        # to a noop (i.e., silently do nothing).
        if entry_point.module_name.split('.') != '__main__':
            util.output_warning(
                'Entry module "{}" basename not "__main__".'.format(
                entry_point.module_name))

        # If this entry module has no entry function, print a non-fatal
        # warning. For unknown reasons, attempting to freeze without an
        # entry function usually causes the frozen executable to reduce to a
        # noop (i.e., silently do nothing).
        if not len(entry_point.attrs):
            util.output_warning(
                'Entry module "{}" entry function undefined.'.format(
                entry_point.module_name))

    # ..................{ SUBCLASS                           }..................
    @abstractmethod
    def _check_frozen_path(self, frozen_pathname: str) -> None:
        '''
        Validate the passed path to which PyInstaller will subsequently write
        the output frozen executable file or directory for the current input
        Python script (in a subclass-specific manner).
        '''
        pass

    @abstractmethod
    def _get_freeze_mode(self) -> str:
        '''
        Get a string constant specific to this subclass.

        This constant should be `file` when freezing in "one-file" mode and
        `dir` when freezing in "one-directory" mode. The BETSE-specific
        `__FREEZE_MODE` environment variable will be set to this constant,
        informing the BETSE-specific PyInstaller specification file run in a
        separate process of which mode to freeze in.
        '''
        pass

    @abstractmethod
    def _get_pyinstaller_options(self) -> list:
        '''
        Get a list of subclass-specific command-line options to be passed to the
        external `pyinstaller` command, when running such command in the absence
        of a previously generated spec file.

        When run with a previously generated spec file, such command effectively
        accepts *no* such options. Hence, this method is only called when no
        such file exists for the current script wrapper to be frozen.
        '''
        pass

class freeze_dir(freeze):
    '''
    Create one platform-specific executable file in one subdirectory of the
    top-level `dist` directory for each previously installed script for the
    current application.

    See Also
    ----------
    freeze
        For further details.
    '''

    description =\
        'freeze all installed scripts to platform-specific executable directories'
    '''
    Command description printed when running `./setup.py --help-commands`.
    '''

    def _check_frozen_path(self, frozen_pathname: str) -> None:
        '''
        Validate that the directory to be generated is *not* an existing file
        (e.g., due to a prior run of the `freeze_file` command).

        Additionally, if such directory exists *and* the user passed option
        `--clean`, such directory will be recursively deleted in a safe manner
        (e.g., *not* following symbolic links outside such directory).
        '''
        util.die_unless_dir_or_not_found(frozen_pathname)

    def _get_freeze_mode(self) -> str:
        return 'dir'

    def _get_pyinstaller_options(self) -> list:
        return [
            '--onedir',
        ]

class freeze_file(freeze):
    '''
    Create one platform-specific executable file in the top-level `dist`
    directory for each previously installed script for the current application.

    See Also
    ----------
    freeze
        For further details.
    '''

    description =\
        'freeze all installed scripts to platform-specific executable files'
    '''
    Command description printed when running `./setup.py --help-commands`.
    '''

    def _check_frozen_path(self, frozen_pathname: str) -> None:
        '''
        Validate that the file to be generated is *not* an existing directory
        (e.g., due to a prior run of the "freeze_dir" command).

        Additionally, if such file exists *and* the user passed option
        `--clean`, such file will be deleted.
        '''
        util.die_unless_file_or_not_found(frozen_pathname)

    def _get_freeze_mode(self) -> str:
        return 'file'

    def _get_pyinstaller_options(self) -> list:
        return [
            '--onefile',
        ]

# --------------------( WASTELANDS                         )--------------------
        # if not is_entry_point:
        #     # Human-readable list of the basenames of such script wrappers.
        #     script_wrapper_basenames = '" and "'.join(
        #         self._metadata['entry_point_basenames'])
        #
        #     # Die, you! Die!
        #     raise DistutilsExecError(
        #         'Entry points "{}" not found. {}'.format(
        #             script_wrapper_basenames, freeze.EXCEPTION_ADVICE))

        #FUXME: Oh, boy. We need to verify here that "betse" has actually been
        #installed under the active Python interpreter. The mere existence of
        #script wrappers is a necessary though insufficient condition! Peruse
        #the main codebase for module detection utility functions.
        #FUXME: Such verification *MUST* be duplicated in the script wrappers as
        #well. We absolutely can't permit this issue to persist in the future.
        #FUXME: Something has gone horribly wrong. To help debug this, let's
        #improve build.py to print the "__file__" attribute of the imported
        #"betse.cli.clicli" module in script wrappers: e.g.,
        #
        #    from betse.cli import clicli
        #
        #    # For debugging purposes, print the absolute path of this module.
        #    print('betse.cli.clicli: ' + clicli.__file__)
        #
        #    # Propagate the return value of such function (hopefully a single-
        #    # byte integer) to the calling shell as this script's exit status.
        #    sys.exit(clicli.main())

    # @abstractmethod
    # def _get_script_spec_basename(self, script_basename: str) -> str:
    #     '''
    #     Get the subclass-specific basename of the PyInstaller spec file
    #     converting the platform-independent Python script with the passed
    #     basename into a platform-specific executable file or directory.
    #
    #     To ensure such spec file is recreated and reused in the same directory
    #     as the top-level setuptools script `setup.py`, such file is returned as
    #     a basename relative to such directory.
    #
    #     This method may also perform subclass-specific validation pertaining to
    #     such basename.
    #     '''
    #     pass

    #FUXME: Refactor this method everywhere to accept no arguments.
    # def _get_script_spec_basename(self, script_basename: str) -> str:
    #     assert isinstance(script_basename, str),\
    #         '"{}" not a string.'.format(script_basename)
    #     return 'dir.spec'
    #     # return '{}.dir.spec'.format(script_basename)
    # def _get_script_spec_basename(self, script_basename: str) -> str:
    #     assert isinstance(script_basename, str),\
    #         '"{}" not a string.'.format(script_basename)
    #     return 'file.spec'
        # return '{}.file.spec'.format(script_basename)

            # self._pyinstaller_spec_filename = path.join('.spec')
            #     pyinstaller_dirname,
            #     self._get_script_spec_basename(script_basename))

#FUXME: We'll almost certainly want to modify the output ".spec" file to
#transparently detect the current OS and modify its behaviour accordingly.
#Happily, ".spec" files appear to be mostly Python. As a trivial example, see:
#    https://github.com/suurjaak/Skyperious/blob/master/packaging/pyinstaller.spec
#Also note the platform-specific instructions at:
#    http://irwinkwan.com/2013/04/29/python-executables-pyinstaller-and-a-48-hour-game-design-compo

                # 'Consider running either the "install" or "symlink" '
                # 'subcommand and trying again. '
                # 'See "README.md" for details.'.format(script_wrapper_basenames)
            #FUXME: The following paths should be shell-quoted. Sadly,
            #we were unable to grok a simple solution, so the current
            #lethargic approach stands.

                # Append all common options.
                # pyinstaller_command.extend(pyinstaller_options_common)

#hasattr(self, 'clean')
    # user_options = freeze.user_options
    # user_options = [
    #     ('clean', 'c',
    #      'clean PyInstaller cache of temporary paths before building'),
    # ]

#  Note that
                # "pyinstaller" supports substantially command-line options
                # under this mode of operation than when passed a Python script.
        # Validate that the file to be generated is *NOT* an existing directory
        # (e.g., due to a prior run of the "freeze_dir" command).
        # util.die_unless_file_or_not_found(path.join('dist', script_basename))
    #FUXME: Call such method above and define below as well.
        # Type of freezing performed by this subclass implementation (e.g.,
        # "file" when performing one-file freezing).
        # freeze_type = self.__class__.__name__[len('freeze_'):]

            # To ensure such file is recreated and reused in the same directory
            # as the top-level "setup.py" script, such file's path is specified
            # as a basename and hence relative to such directory.
            # script_spec_basename = '{}.{}.spec'.format(
            #     script_basename, freeze_type)
#FUXME: Also make a "freeze_dir" class. Since such class is intended only for
#debugging, such class' run() method should also pass the "--debug" option to
#"pyinstaller".

        # Copy the "install_dir" attribute from the existing "install_scripts"
        # attribute of a temporarily instantiated "symlink" object.
        #
        # Why? Because setuptools.
    #FUXME: Does this actually work? If so, replicate to the other modules in
    #this package as well.

            # Basename of the PyInstaller ".spec" file converting such
            # platform-independent script into a platform-specific executable.
            # Since such file is specific to both the basename of such script
            # *AND* the type of the current operating system, such substrings
            # are embedded in such filename.
            # script_spec_basename = '{}.{}.spec'.format(
            #     script_basename, util.get_os_type())

                # # Basename of the currently generated spec file.
                # script_spec_basename_current = '{}.spec'.format(script_basename)
                #
                # # Rename such file to have the operating system-specific
                # # basename expected by the prior conditional (on the next
                # # invocation of this setuptools command).
                # #
                # # Note that "pyinstaller" accepts an option "--name" permitting
                # # the basename of such file to be specified prior to generating
                # # such file. Unfortunately, such option *ALSO* specifies the
                # # basename of the generated executable. Since we only
                # util.move_file(
                #     script_spec_basename_current, script_spec_basename)

                #FUXME: If a ".spec" file exists, such file should be passed rather
                #than such script's basename. Note that, when passing such file,
                #most CLI options are ignored; of those we use below, only
                #"--noconfirm" appears to still be respected.

                    # path.basename(script_spec_basename))
            # Absolute path of the PyInstaller ".spec" file converting such
            # platform-independent script into a platform-specific executable.
            # Since such file is specific to both the basename of such script
            # *AND* the type of the current operating system, such substrings
            # are embedded in such filename.
            # self._pyinstaller_spec_filename = path.join(
            #     script_basename, util.get_os_type()
            # )

            # util.output_sans_newline(
            #     'Searching for existing spec file "{}"... '.format(
            #         path.basename(self._pyinstaller_spec_filename)))
#".app"-suffixed
        # # List of shell words common to all "pyinstaller" commands called below.
        # command_words_base = [
        #     'pyinstaller',
        #     '--onefile',
        #     # Overwrite existing output paths under the "dist/" subdirectory
        #     # without confirmation, the default behaviour.
        #     '--noconfirm',
        # ]
            # command_words_base
            #
            # if script_type == 'console':
            #     --console
                # 'Disabling compression of output executables.'
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
