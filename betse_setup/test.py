#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Application-specific ``test`` subcommand for :mod:`setuptools`.
'''

# ....................{ IMPORTS                            }....................
from betse.util.type.obj import objects
from betse_setup import buputil
from distutils.errors import DistutilsClassError
from setuptools import Command

# ....................{ COMMANDS                           }....................
def add_setup_commands(metadata: dict, setup_options: dict) -> None:
    '''
    Add the ``test`` subcommand to the passed dictionary of :mod:`setuptools`
    options.
    '''

    buputil.add_setup_command_classes(metadata, setup_options, test)

# ....................{ CLASSES ~ base                     }....................
class test(Command):
    '''
    Command class testing the current application with :mod:`pytest`.

    Attributes
    ----------
    match_name : str
        Python-evaluatable expression (e.g., ``'test_method or test other'``)
        conditionally matching a substring of the name of each test function or
        class to be run if any *or* `None` if all tests are to be run
        unconditionally.
    _pytest_public : PackageType
        Top-level public :mod:`pytest` package, imported by the :meth:`run`
        method.
    _pytest_private : PackageType
        Top-level private :mod:`_pytest` package, imported by the :meth:`run`
        method.

    See Also
    ----------
    https://github.com/pytest-dev/pytest-runner
        Recommended implementation of a :mod:`setuptools` command invoking
        :mod:`pytest`, which is overly obsessed with egg-based local
        installation of test dependencies rather than actual testing.
        Naturally, this is less than useful.
    https://pytest.org/latest/goodpractices.html#manual-integration
        Fallback implementation for a :mod:`setuptools` command invoking
        :mod:`pytest`, which is sufficiently trivial as to be effectively
        useless. This is also less than useful.
    '''

    # ..................{ ATTRIBUTES                         }..................
    description = 'run py.test-driven functional and unit tests'
    '''
    Command description printed on running ``./setup.py --help-commands``.
    '''


    user_options = [
        # Boolean options.
        (
            'no-capture', 's',
            'Prevent py.test from silently capturing any '
            'standard error, standard output, or logging messages '
            'emitted by tests. '
            '(By default, py.test silently captures all three.) '
            'All three will be written in real-time "as is" to their '
            'respective file handles (e.g., current terminal, logfile).'
        ),

        # Argument options.
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
        (
            'export-sim-conf-dir=', None,
            'Target directory into which all '
            'source simulation configuration directories produced by '
            '"@skip_unless_export_sim_conf"-marked tests are to be copied.'
        ),
    ]
    '''
    List specifying command-line options accepted by this command.

    Each item of this list is a 3-tuple ``(long_name, short_name, synopsis)``
    specifying a single option, where:

    * ``long_name`` is this option's mandatory GNU-style long form excluding
      ``--`` prefix. Additionally:
      * If this option accepts a mandatory argument, this long form should be
        suffixed by ``=``.
      * Else, this long form should *not* be suffixed by ``=``.
    * ``short_name`` is either:
      * This option's option *NIX-style short form excluding ``-`` prefix.
      * ``None`` if this option has no such form.
    * ``synopsis`` is this option's mandatory help documentation.

    For each such option, an attribute of the same name as this option's long
    form *must* be explicitly initialized in the :meth:`initialize_options`
    method. :mod:`setuptools` fails to recognize options for which this is
    *not* the case. (You fail a simple sanity check yet again,
    :mod:`setuptools`.)

    See Also
    ----------
    http://ilostmynotes.blogspot.ca/2009/04/python-distutils-installer-and.html
        Inarguably, the best (albeit unofficial) documentation on this list.
    '''

    # ..................{ SUPERCLASS                         }..................
    def initialize_options(self):
        '''
        Declare option-specific attributes subsequently initialized by
        :meth:`finalize_options`.

        If this function is *not* defined, the default implementation of this
        method raises an inscrutable :mod:`distutils` exception. If these
        attributes are *not* declared, the subsequent call to
        :meth:`self.set_undefined_options` raises an inscrutable
        :mod:`setuptools` exception. (This is terrible. So much hate.)
        '''

        # Option-specific public attributes. For each option declared by the
        # "user_options" list above, a public attribute of the same name as this
        # option's long form *MUST* be initialized here to its default value.
        self.no_capture = None
        self.match_name = None
        self.export_sim_conf_dir = None

        # setuptools-specific public attributes.
        # self.install_dir = None

        # Custom private attributes.
        self._pytest_public = None
        self._pytest_private = None


    def finalize_options(self):
        '''
        Default undefined command-specific options to the options passed to the
        current parent command if any (e.g., ``symlink``).
        '''

        pass


    def run(self):
        '''
        Run the current command and all subcommands thereof.
        '''

        # Import and monkey-patch "py.test" *BEFORE* running "py.test".
        self._init_pytest()
        self._patch_pytest()
        self._run_pytest()

    # ..................{ PRIVATE                            }..................
    def _init_pytest(self) -> None:
        '''
        Dynamically import the top-level :mod:`pytest` package into this
        object's :attr:`_pytest` instance variable.

        Specifically, this method imports the top-level `_pytest.main` module,
        providing programmatic access to py.test's CLI implementation. While
        py.test is also runnable as an external command, doing so invites
        non-trivial complications. Unlike most Python applications (e.g.,
        PyInstaller), py.test is *not* dynamically importable via the following
        import machinery:

            # Don't do this.
            >>> pytest_main = buputil.import_module(
            ...    'pytest.main', exception_message=(
            ...    'py.test not installed under the current Python interpreter.'))

        Attempting to do so reliably raises exceptions resembling:

            AttributeError: 'module' object has no attribute '__path__'

        The reason why appears to be package ``__path__`` manipulations
        dynamically performed by the private :mod:`_pytest` package. Sadly,
        attempting to import :func:`_pytest.main` fails with an even more
        inscrutable exception. The solution, of course, is to brute-force it
        by only dynamically importing the root :mod:`pytest` package.
        '''

        # Import the public "pytest" package *BEFORE* the private "_pytest"
        # package. While importation order is typically ignorable, imports can
        # technically have side effects. Tragicomically, this is the case here.
        # Importing the public "pytest" package establishes runtime
        # configuration required by submodules of the private "_pytest" package.
        # The former *MUST* always be imported before the latter. Failing to do
        # so raises obtuse exceptions at runtime... which is bad.
        self._pytest_public = buputil.import_module(
            'pytest', exception_message=(
            'py.test not installed under the current Python interpreter.'))
        self._pytest_private = buputil.import_module(
            '_pytest', exception_message=(
            'py.test not installed under the current Python interpreter.'))


    def _patch_pytest(self) -> None:
        '''
        Monkey-patch the :mod:`pytest` framework in the active Python
        interpreter, altering startup :mod:`pytest` functionality in an
        early-time manner *not* permitted within :mod:`pytest` plugins.

        `pytest` plugins (e.g., ``conftest`` submodules of a test suite) are
        imported by :mod:`pytest` *after* :mod:`pytest` startup and hence cannot
        be used to alter startup :mod:`pytest` functionality in an early-time
        manner.

        Specifically, this method monkey-patches:

        * The :meth:`CaptureManager._getcapture` method to capture stderr but
          *not* stdout (rather than neither stderr nor stdout) when the
          ``py.test`` command is passed either the ``-s`` or ``--capture=no``
          CLI options. The default approach of *not* capturing stderr prevents
          ``py.test`` from capturing and hence reporting error messages in
          failure reports, requiring tedious upwards scrolling through test
          output to find the corresponding error messages.
        '''

        # py.test classes to be monkey-patched.
        from _pytest.capture import CaptureManager, FDCapture, MultiCapture

        # If the private method to be monkey-patched no longer exists, py.test
        # is either broken or unsupported. In either case, raise an exception.
        if not objects.is_method(CaptureManager, '_getcapture'):
            raise DistutilsClassError(
                'Method pytest.capture.CaptureManager._getcapture() '
                'not found. The current version of py.test is either '
                'broken (unlikely) or unsupported (likely).'
            )

        # Old method to be monkey-patched.
        _getcapture_old = CaptureManager._getcapture

        # New method applying this monkey-patch.
        def _getcapture_new(self, method):
            if method == "no":
                return MultiCapture(
                    out=False, err=True, in_=False, Capture=FDCapture)
            else:
                return _getcapture_old(self, method)

        # Replace the old with the new method.
        CaptureManager._getcapture = _getcapture_new


    def _run_pytest(self) -> None:
        '''
        Call the :func:`pytest.main` function in the active Python interpreter,
        passed CLI options corresponding to the CLI options passed to this
        setuptools command.
        '''

        # List of all shell words to be passed as arguments to "py.test".
        #
        # Whereas the top-level "pytest.ini" configuration file lists options
        # unconditionally passed to *ALL* "py.test" invocations, this list only
        # lists options conditionally applicable to invocations run within the
        # "python3 setup.py test" subcommand.
        #
        # For example, the "--capture=no" option should be conditionally enabled
        # *ONLY* when applying the monkey-patch applied by this setuptools
        # subcommand (namely, here). This option is intentionally omitted from
        # the top-level "pytest.ini" configuration file.
        pytest_args = [
            # When testing interactively, prevent py.test from capturing stdout
            # but *NOT* stderr. By default, py.test captures and delays printing
            # stdout until after test completion. While a possibly suitable
            # default for short-lived unit tests, such capturing is unsuitable
            # for long-lived functional tests.
            #
            # Note that this option is monkey-patched by the _patch_pytest()
            # method to capture only stdout. By default, this option captures
            # neither stdout (which is good) nor stderr (which is bad).
            '--capture=no',

            # When testing interactively, halt testing on the first failure.
            # Permitting multiple failures complicates failure output,
            # especially when every failure after the first is a result of the
            # same underlying issue.
            #
            # When testing non-interactively, testing is typically *NOT* halted
            # on the first failure. Hence, this option is confined to this
            # subcommand rather than added to the general-purpose "pytest.ini"
            # configuration.
            '--maxfail=1',
        ]

        #FIXME: Disable "xdist" if at least one serial test exists. Currently,
        #none do. Theoretically, they could. They once did and certaily could
        #again. Serial functional tests assume that all tests to be serialized
        #are run in the same Python process. "pytest-xdist" currently provides
        #no means of doing so; instead, "pytest-xdist" assigns all tests to
        #arbitrary test slaves and hence Python processes. Until "pytest-xdist"
        #permits tests to be isolated to the same test slave, "pytest-xdist"
        #must *NOT* be enabled. Non-fatal warning to output might resemble:
        #    buputil.output_warning(
        #        'py.test plugin "pytest-xdist" fails to support serial tests.')
        #    buputil.output_warning(
        #        'Tests will *NOT* be parallelized across multiple processors.')

        # True only if the optional third-party "pytest-xdist" plugin is both
        # importable and *NOT* explicitly disabled below (e.g., due to the end
        # user having passed CLI options incompatible with this plugin).
        is_xdist = buputil.is_module('xdist')

        #FIXME: Disabled for the moment. "xdist" appears to be inexplicably
        #failing with non-human-readable exceptions. The lack of official
        #support for test parallelization in py.test is becoming clear.
        #
        #When "xdist" is (eventually) stabilized, excise this reassignment here
        #*AND* reenable the conditional below.
        is_xdist = False

        #FIXME: Reenable this conditional.
        # If this plugin is unimportable, print a non-fatal warning. Due to the
        # cost of running tests, parallelization is highly recommended.
        # if not is_xdist:
        #     buputil.output_warning(
        #         'Optional py.test plugin "pytest-xdist" not found.')
        #     buputil.output_warning(
        #         'Tests will *NOT* be parallelized across multiple processors.')

        # Pass options passed to this subcommand to this py.test command,
        # converting long option names specific to this subcommand (e.g.,
        # "--no-capture") to short option names recognized by py.test (e.g.,
        # "-s"). Sadly, py.test typically recognizes only the latter.
        #
        # If the "-s" option is passed...
        if self.no_capture is not None:
            # If "pytest-xdist" is importable, print a non-fatal warning.
            # This plugin silently ignores all capture redirection CLI options,
            # including "-s", "--no-capture", and "--capture" options. It
            # appears unlikely that this will ever be solved. For details, see:
            # https://github.com/pytest-dev/pytest/issues/680
            if is_xdist:
                # Print a non-fatal warning.
                buputil.output_warning(
                    'Option "-s" unsupported by py.test plugin "pytest-xdist".')
                buputil.output_warning(
                    'Tests will *NOT* be parallelized across '
                    'multiple processors.')

                # Disable this plugin. While parallelization is important,
                # respecting end user wishes with respect to capture redirection
                # is paramount and hence takes precedence.
                is_xdist = False

            # Pass the "-s" option to py.test.
            pytest_args.append('-s')
            # pytest_args.append('--capture=no')

        # If any remaining option is passed, forward these option unmodified
        # onto "py.test". For safety, these option's arguments are intentionally
        # *NOT* shell-quoted. Doing so unnecessarily adds an additional level of
        # quoting... which is bad.
        if self.match_name is not None:
            pytest_args.extend(('-k', self.match_name,))
        if self.export_sim_conf_dir is not None:
            pytest_args.extend(
                ('--export-sim-conf-dir', self.export_sim_conf_dir,))

        # If "pytest-xdist" is both importable and *NOT* explicitly disabled...
        if is_xdist:
            # Notify the user of test parallelization.
            print('Optional py.test plugin "pytest-xdist" found.')
            print('Tests will be parallelized across all available processors.')

            # Instruct "pytest-xdist" to autodetect and parallelize tests to all
            # available processors.
            pytest_args.extend(['-n', 'auto'])

        # Instruct "py.test" of the relative pathname to the top-level directory
        # for this project. On startup, "py.test" internally:
        #
        # * Sets its "rootdir" property to this pathname in absolute form.
        # * Sets its "inifile" property to the concatenation of this pathname
        #   with the basename "pytest.ini" if this top-level configuration file
        #   exists.
        # * Prints the initial values of these properties to stdout.
        #
        # *THIS IS ESSENTIAL.* If "py.test" is *NOT* explicitly passed this
        # relative pathname as an argument, "py.test" typically fails to set
        # these properties to the expected paths. For unknown reasons
        # (presumably unresolved "py.test" issues), "py.test" instead sets
        # "rootdir" to the absolute path of the current user's home directory
        # and "inifile" to "None" -- because no current user's home directory
        # contains a "pytest.ini" file. The error output resembles:
        #
        #    $ ./test -k test_sim_export --export-sim-conf-dir ~/tmp/yolo
        #    running test
        #    Running py.test with arguments: ['--capture=no', '--maxfail=1', '-k', 'test_sim_export', '--export-sim-conf-dir', '/home/leycec/tmp/yolo']
        #    usage: setup.py [options] [file_or_dir] [file_or_dir] [...]
        #    setup.py: error: unrecognized arguments: --export-sim-conf-dir
        #      inifile: None
        #      rootdir: /home/leycec
        #
        # See the following official documentation for further details, entitled
        # "Initialization: determining rootdir and inifile":
        #     https://docs.pytest.org/en/latest/customize.html
        pytest_args.append('.')

        # Fork the "py.test" command, propagating its exit status as our own.
        print('Running py.test with arguments: {}'.format(pytest_args))
        buputil.exit_with_status(self._pytest_public.main(pytest_args))
