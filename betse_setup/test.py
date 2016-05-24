#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
BETSE-specific `test` subcommand for `setuptools`.
'''

# ....................{ IMPORTS                            }....................
from betse_setup import util
from setuptools import Command

# ....................{ COMMANDS                           }....................
def add_setup_commands(metadata: dict, setup_options: dict) -> None:
    '''
    Add the `test` subcommand to the passed dictionary of `setuptools` options.
    '''
    util.add_setup_command_classes(metadata, setup_options, test)

# ....................{ CLASSES ~ base                     }....................
class test(Command):
    '''
    Command class testing the current application with `py.test`.

    Attributes
    ----------
    match_name : str
        Python-evaluatable expression (e.g., `'test_method or test other'`)
        conditionally matching a substring of the name of each test function or
        class to be run if any _or_ `None` if all tests are to be run
        unconditionally.
    _pyinstaller_command : list
        List of all shell words of the PyInstaller command to be run.

    See Also
    ----------
    https://github.com/pytest-dev/pytest-runner
        Recommended implementation of a `setuptools` command invoking `pytest`,
        which is overly obsessed with egg-based local installation of test
        dependencies rather than actual testing. This is less than useful.
    https://pytest.org/latest/goodpractices.html#manual-integration
        Fallback implementation for a `setuptools` command invoking `pytest`,
        which is sufficiently trivial as to be effectively useless.
    '''

    description = 'run py.test-driven functional and unit tests'
    '''
    Command description printed on running `./setup.py --help-commands`.
    '''

    user_options = [
        (
            'no-capture', 's',
            'Prevent py.test from silently capturing any '
            'standard error, standard output, or logging messages '
            'emitted by tests. '
            '(By default, py.test silently captures all three.) '
            'All three will be written in real-time "as is" to their '
            'respective file handles (e.g., current terminal, logfile).'
        ),
        (
            'match-name=', 'k',
            'Only run tests which match the given substring expression. '
            'An expression is a python evaluatable expression '
            'where all names are substring-matched against '
            'test names and their parent classes. '
            'Example: -k "test_method or test other" matches '
            'all test functions and classes whose names contain '
            'either "test_method" or "test_other". '
            'Keywords are also matched to classes and functions containing '
            'extra names in their "extra_keyword_matches" set as well as to '
            'functions which have names assigned directly to them.'
        ),
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
        Inarguably, the best (albeit unofficial) documentation on this list.
    '''

    # ..................{ SUPERCLASS                         }..................
    def initialize_options(self):
        '''
        Declare option-specific attributes subsequently initialized by
        `finalize_options()`.

        If this function is _not_ defined, the default implementation of this
        method raises an inscrutable `distutils` exception. If these attributes
        are _not_ declared, the subsequent call to
        `self.set_undefined_options()` raises an inscrutable `setuptools`
        exception. (This is terrible. So much hate.)
        '''

        # Option-specific public attributes. For each option declared by the
        # "user_options" list above, a public attribute of the same name as this
        # option's long form *MUST* be initialized here to its default value.
        self.no_capture = None
        self.match_name = None

        # setuptools-specific public attributes.
        # self.install_dir = None

        # Custom private attributes.
        # self._pyinstaller_command = None


    def finalize_options(self):
        '''
        Default undefined command-specific options to the options passed to the
        current parent command if any (e.g., `symlink`).
        '''
        pass


    def run(self):
        '''Run the current command and all subcommands thereof.'''

        # List of all shell words to be passed as arguments to py.test.
        pytest_args = []

        # Pass options passed to this subcommand to this py.test command,
        # converting long option names specific to this subcommand (e.g.,
        # "--no-capture") to short option names recognized by py.test (e.g.,
        # "-s"). Sadly, py.test typically recognizes only the latter.
        if self.no_capture is not None:
            pytest_args.append('-s')
        if self.match_name is not None:
            pytest_args.extend(['-k', util.shell_quote(self.match_name)])

        # If the optional third-party "pytest-xdist" plugin is installed, pass
        # options specific to this plugin implicitly parallelizing tests to a
        # number of processors (hopefully) autodetected at runtime.
        if util.is_module('xdist'):
            pytest_args.extend(['-n', 'auto'])
        # Else, print a non-fatal warning. Due to the cost of running tests,
        # parallelization is highly recommended.
        else:
            util.output_warning(
                'Optional py.test plugin "pytest-xdist" not found.')
            util.output_warning(
                'Tests will *NOT* be parallelized across multiple processors.')

        # py.test's top-level "main" module, providing programmatic access to
        # its CLI implementation. While py.test is also runnable as an external
        # command, doing so invites non-trivial complications. Unlike most
        # Python applications (e.g., PyInstaller), py.test is *NOT* dynamically
        # importable via the following import machinery:
        #
        #     # Don't do this.
        #     pytest_main = util.import_module(
        #         'pytest.main', exception_message=(
        #         'py.test not installed under the current Python interpreter.'))
        #
        # Attempting to do so reliably raises exceptions resembling:
        #
        #     AttributeError: 'module' object has no attribute '__path__'
        #
        # The reason why appears to be package "__path__" manipulations
        # dynamically performed by the private "_pytest" package. Sadly,
        # attempting to import "_pytest.main" fails with an even more
        # inscrutable exception. The solution, of course, is to brute-force it
        # by only dynamically importing the root "pytest" package. *SHRUG*
        pytest = util.import_module(
            'pytest', exception_message=(
            'py.test not installed under the current Python interpreter.'))

        # Run py.test, propagating its exit status as our own up to the caller.
        print('Running py.test with arguments: {}'.format(pytest_args))
        util.exit_with_status(pytest.main(pytest_args))
