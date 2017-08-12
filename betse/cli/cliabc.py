#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Abstract base classes for defining command line interface (CLI) applications.
'''

# ....................{ IMPORTS                            }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To raise human-readable exceptions on application startup, the
# top-level of this module may import *ONLY* from submodules guaranteed to:
# * Exist, including standard Python and application modules.
# * Never raise exceptions on importation (e.g., due to module-level logic).
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

import sys
from abc import ABCMeta, abstractmethod
from betse import ignition, pathtree
from betse.cli import cliutil
from betse.cli.cliopt import (
    CLIOptionArgEnum,
    CLIOptionArgStr,
    CLIOptionBoolTrue,
    CLIOptionVersion,
)
from betse.lib import libs
from betse.util.io.log import logs, logconfig
from betse.util.io.log.logenum import LogLevel
from betse.util.path.command import cmds
from betse.util.path.command.cmdarg import HelpFormatterParagraph
from betse.util.path.command.cmdexit import SUCCESS, FAILURE_DEFAULT
from betse.util.py.pyprof import profile_callable, ProfileType
from betse.util.type import types
from betse.util.type.types import (
    type_check, ArgParserType, MappingType, SequenceOrNoneTypes)

# ....................{ SUPERCLASS                         }....................
class CLIABC(object, metaclass=ABCMeta):
    '''
    Abstract base class of all top-level command line interface (CLI)
    subclasses, suitable for use by both CLI and GUI front-ends for BETSE.

    Attributes
    ----------
    _arg_list : list
        List of all passed command-line arguments as unparsed raw strings.
    _arg_parser : ArgParserType
        `argparse`-specific parser of command-line arguments.
    _arg_parser_kwargs : dict
        Dictionary of keyword arguments which which to create argument parsers,
        suitable for passing to both the :meth:`ArgParserType.__init__` and
        :meth:`ArgParserType.add_parser` methods. Since the :mod:`argparse` API
        provides multiple methods rather than a single method for creating
        argument parsers, this versatile dictionary is preferred over a
        monolithic factory-based approach (e.g., a `_make_arg_parser()` method).
    _args : argparse.Namespace
        :mod:`argparse`-specific container of all passed command-line arguments.
        See "Attributes (_args)" below for further details.
    _profile_filename : str
        Absolute or relative path of the dumpfile to export a profile of the
        current execution to if :attr:`_profile_type` is
        :attr:`ProfileType.CALL` *or* ignored otherwise.
    _profile_type : ProfileType
        Type of profiling to be performed if any.

    Attributes (of :attr:`_args`)
    ----------
    is_verbose : bool
        ``True`` only if low-level debugging messages are to be logged. Defaults
        to ``False``.
    log_filename : str
        Absolute or relative path of the file to log to. Defaults to the
        absolute path of BETSE's default user-specific logfile.
    log_level : str
        Minimum level of messages to be logged to :attr:`log_filename`,
        formatted as the lowercased name of a :class:`LogLevel` enumeration
        member. Defaults to ``info``.
    stderr_level : str
        Minimum level of messages to be redirected to stderr, formatted as the
        lowercased name of a :class:`LogLevel` enumeration member. Defaults to
        ``warning``.
    stdout_level : str
        Minimum level of messages to be redirected to stdout, formatted as the
        lowercased name of a :class:`LogLevel` enumeration member. Defaults to
        ``info``.
    profile_filename : str
        Absolute or relative path of the dumpfile to export a profile of the
        current execution to if :attr:`profile_type` is ``call`` *or* ignored
        otherwise.  Defaults to the absolute path of BETSE's default
        user-specific profile dumpfile.
    profile_type : str
        Type of profiling to be performed if any, formatted as the lowercased
        name of a :class:`ProfileType` enumeration member. Defaults to ``none``.
    '''

    # ..................{ INITIALIZERS                       }..................
    def __init__(self):

        # Initialize subclasses performing diamond inheritance if any.
        super().__init__()

        # For safety, nullify all remaining attributes.
        self._arg_list = None
        self._arg_parser = None
        self._arg_parser_kwargs = None
        self._args = None
        self._profile_filename = None
        self._profile_type = None

    # ..................{ RUNNERS                            }..................
    @type_check
    def run(self, arg_list: SequenceOrNoneTypes = None) -> int:
        '''
        Run the command-line interface (CLI) defined by this subclass with the
        passed argument list if non-``None`` *or* the external argument list
        passed on the command line (i.e., :data:`sys.argv`) otherwise.

        Parameters
        ----------
        arg_list : optional[SequenceTypes]
            Sequence of zero or more arguments to pass to this interface.
            Defaults to ``None``, in which case arguments passed on the command
            line (i.e., :data:`sys.argv`) are used instead.

        Returns
        ----------
        int
            Exit status of this interface in the range ``[0, 255]``.
        '''

        # Default unpassed arguments to those passed on the command line,
        # ignoring the first element of "sys.argv" (i.e., the filename of the
        # command from which the current Python process was spawned).
        if arg_list is None:
            # logs.log_info('Defaulting to sys.argv')
            arg_list = sys.argv[1:]
        assert types.is_sequence_nonstr(arg_list), (
            types.assert_not_sequence_nonstr(arg_list))
        # print('BETSE arg list (in run): {}'.format(arg_list))

        # Classify arguments for subsequent use.
        self._arg_list = arg_list

        try:
            # (Re-)initialize this application *BEFORE* performing subsequent
            # logic assuming this application to have already been initialized.
            self._ignite_app()

            # Parse these arguments *AFTER* initializing logging, ensuring
            # logging of exceptions raised by this parsing.
            self._parse_args()

            # (Re-)initialize all mandatory runtime dependencies *AFTER* parsing
            # and handling all logging-specific CLI options and hence finalizing
            # the logging configuration for the active Python process. This
            # initialization integrates the custom logging and debugging schemes
            # implemented by these dependencies with that implemented by BETSE.
            libs.reinit(
                matplotlib_backend_name=self._args.matplotlib_backend_name)

            # Run the command-line interface (CLI) defined by this subclass,
            # profiled by the type specified by the "--profile-type" option.
            profile_callable(
                call=self._do,
                profile_type=self._profile_type,
                profile_filename=self._profile_filename,
            )
            # libs.die_unless_runtime_optional('networkx')
            # raise ValueError('Test exception handling.')

            # Exit with successful exit status from the current process.
            return SUCCESS
        except Exception as exception:
            # Handle this exception.
            self._handle_exception(exception)

            # Exit with failure exit status from the current process. If this
            # exception provides a system-specific exit status, use this status;
            # else, use the default failure status (i.e., 1).
            #
            # Ignore the Windows-specific "winerror" attribute provided by
            # "WindowsError"-based exceptions. While more fine-grained than the
            # "errno" attribute, "winerror" values are *ONLY* intended to be
            # used internally rather than returned as an exit status.
            return getattr(exception, 'errno', FAILURE_DEFAULT)

    # ..................{ ARGS                               }..................
    def _parse_args(self) -> None:
        '''
        Parse all currently passed command-line arguments.

        In order, this method:

        * Creates and configures an argument parser with sensible defaults.
        * Calls the subclass-specific :meth:`_config_arg_parsing` method,
          defaulting to a noop.
        * Parses all arguments with this parser.
        '''

        # Create and configure all argument parsers.
        self._init_arg_parsers()

        # Parse unnamed string arguments into named, typed arguments.
        self._args = self._arg_parser.parse_args(self._arg_list)

        # Parse top-level options globally applicable to *ALL* subcommands.
        self._parse_options_top()


    def _init_arg_parsers(self) -> None:
        '''
        Create and configure all argument parsers, including both the top-level
        argument parser and all subparsers of that parser.
        '''

        # Create and classify the top-level argument parser.
        self._init_arg_parser_top()

        # Configure subclass-specific argument parsing.
        self._config_arg_parsing()


    def _init_arg_parser_top(self) -> None:
        '''
        Create and classify the top-level argument parser.
        '''

        # Initialize these keyword arguments.
        self._arg_parser_kwargs = {
            # Wrap non-indented lines in help and description text as paragraphs
            # while preserving indented lines in such text as is.
            'formatter_class': HelpFormatterParagraph,
        }

        # Dictionary of keyword arguments initializing the core argument parser.
        arg_parser_top_kwargs = {
            # Script name.
            'prog': cmds.get_current_basename(),
        }

        # Update this dictionary with preinitialized arguments.
        arg_parser_top_kwargs.update(self._arg_parser_kwargs)

        # Update this dictionary with subclass-specific arguments.
        arg_parser_top_kwargs.update(self._arg_parser_top_kwargs)

        # Core argument parser.
        self._arg_parser = ArgParserType(**arg_parser_top_kwargs)

        # For each top-level option, add an argument parsing this option to this
        # argument subparser.
        for option in self._make_options_top():
            option.add(self._arg_parser)

    # ..................{ ARGS ~ options                     }..................
    def _make_options_top(self) -> tuple:
        '''
        Tuple of all :class:`CLIOptionABC` instances defining the top-level
        CLI options accepted by this application.

        For each such option, the :meth:`_init_arg_parser_top` method adds a
        corresponding argument to the top-level argument parser (i.e.,
        :attr:`_arg_parser`).

        Design
        ----------
        Subclasses requiring subclass-specific options are encouraged to:

        * Override this method by:
          * Calling the superclass implementation.
          * Extending the returned tuple with all subclass-specific options.
        * Override the :meth:`_parse_options_top` method by:
          * Calling the superclass implementation.
          * Handling all subclass-specific instance variables parsed into the
            :attr:`self._args` container from these options.

        Caveats
        ----------
        Order is significant, defining the order that the ``betse --help``
        command synopsizes these options in. Options *not* listed here are
        *not* parsed by argument subparsers and hence effectively ignored.

        Returns
        ----------
        tuple
            Tuple of all such :class:`CLIOptionABC` instances.
        '''

        # Singleton logging configuration for the current Python process.
        log_config = logconfig.get()

        # Return a tuple of all default top-level options.
        return (
            CLIOptionBoolTrue(
                short_name='-v',
                long_name='--verbose',
                synopsis='print and log all messages verbosely',
            ),

            CLIOptionVersion(
                short_name='-V',
                long_name='--version',
                synopsis='print program version and exit',
                version=cliutil.get_version(),
            ),

            CLIOptionArgStr(
                long_name='--matplotlib-backend',
                synopsis=(
                    'name of matplotlib backend to use '
                    '(see: "betse info")'
                ),
                var_name='matplotlib_backend_name',
                default_value=None,
            ),

            CLIOptionArgStr(
                long_name='--log-file',
                synopsis=(
                    'file to log to '
                    '(defaults to "{default}") '
                ),
                var_name='log_filename',
                default_value=log_config.filename,
            ),

            CLIOptionArgEnum(
                long_name='--log-level',
                synopsis=(
                    'minimum level of messages to log to "--log-file" '
                    '(defaults to "{default}") '
                    '[overridden by "--verbose"]'
                ),
                enum_type=LogLevel,
                enum_default=log_config.file_level,
            ),

            CLIOptionArgEnum(
                long_name='--profile-type',
                synopsis='''
    type of profiling to perform (defaults to "{default}"):
    ;* "none", disabling profiling
    ;* "call", profiling callables (functions, methods)
    ;* "size", profiling object sizes (requires "pympler")
    ''',
    # ;* "line", profiling code lines (requires "pprofile")
                enum_type=ProfileType,
                enum_default=ProfileType.NONE,
            ),

            CLIOptionArgStr(
                long_name='--profile-file',
                synopsis=(
                    'file to profile to unless "--profile-type=none" '
                    '(defaults to "{default}")'
                ),
                var_name='profile_filename',
                default_value=pathtree.get_profile_default_filename(),
            ),
        )


    def _parse_options_top(self) -> None:
        '''
        Parse top-level options globally applicable to all subcommands.

        Design
        ----------
        Subclasses requiring subclass-specific options are encouraged to
        override this method. See the :meth:`_make_options_top` method for
        further details.
        '''

        # Configure logging options *BEFORE* all remaining options, ensuring
        # proper logging of messages emitted by the latter.
        self._parse_options_top_log()
        self._parse_options_top_profile()


    def _parse_options_top_log(self) -> None:
        '''
        Parse top-level logging options globally applicable to all subcommands.
        '''

        # Singleton logging configuration for the current Python process.
        log_config = logconfig.get()

        # Configure logging according to the passed options. Note that order of
        # assignment is insignificant here.
        # print('is verbose? {}'.format(self._args.is_verbose))
        log_config.is_verbose = self._args.is_verbose
        log_config.filename = self._args.log_filename
        log_config.file_level = LogLevel[self._args.log_level.upper()]

        # Display a human-readable synopsis of this application.
        self._show_header()

        # Log all string arguments passed to this command.
        logs.log_debug('Passed argument list: {}'.format(self._arg_list))


    def _parse_options_top_profile(self) -> None:
        '''
        Parse top-level profiling options globally applicable to all
        subcommands.
        '''

        # Classify the passed profiling options, converting the profiling type
        # from a raw lowercase string into a full-fledged enumeration member.
        # See _parse_options_top_log() for further detail on this conversion.
        self._profile_filename = self._args.profile_filename
        self._profile_type = ProfileType[self._args.profile_type.upper()]

    # ..................{ IGNITERS                           }..................
    def _ignite_app(self) -> None:
        '''
        (Re-)initialize this application *BEFORE* performing subsequent logic
        assuming this application to have already been initialized.

        Defaults to (re-)initializing all low-level BETSE logic. Subclasses may
        override this method to perform additional initialization, in which
        case this superclass method should still be called to initialize BETSE.
        '''

        # (Re-)initialize BETSE. Note that calling the ignition.init() function:
        #
        # * Suffices when BETSE is *NOT* running under a test suite.
        # * Fails to suffice if BETSE is running under a test suite, in
        #   which case this suite may run each test from within the same
        #   Python process. Due to caching internally performed by the
        #   ignition.init() function, calling that function here would fail
        #   to re-initialize BETSE in any test except the first. To
        #   model the real world as closely as reasonable, the
        #   ignition.reinit() function is called instead.
        ignition.reinit()

    # ..................{ EXCEPTIONS                         }..................
    @type_check
    def _handle_exception(self, exception: Exception) -> None:
        '''
        Handle the passed uncaught exception, typically by at least logging and
        optionally displaying this exception's detailed traceback to end users.

        Defaults to merely logging this exception. Subclasses may override this
        method to perform additional exception handling, in which case this
        superclass method should still be called to log exceptions.
        '''

        # Log this exception.
        logs.log_exception(exception)

    # ..................{ SUBCLASS ~ mandatory               }..................
    # The following methods *MUST* be implemented by subclasses.

    @abstractmethod
    def _do(self) -> object:
        '''
        Implement this command-line interface (CLI) in a subclass-specific
        manner, returning an arbitrary object produced by this logic to be
        memory profiled when the ``--profile-type=size`` CLI option is passed.
        '''

        pass

    # ..................{ SUBCLASS ~ optional                }..................
    # The following methods may but need *NOT* be implemented by subclasses.

    @property
    def _arg_parser_top_kwargs(self) -> MappingType:
        '''
        Subclass-specific dictionary of all keyword arguments to be passed to
        the :meth:`ArgP"arserType.__init__` method of the top-level argument
        parser for this CLI.

        Defaults to the empty dictionary.

        See Also
        ----------
        :meth:`_init_arg_parser_top`
            Method initializing this argument parser with this dictionary *and*
            keyword arguments common to all argument parsers.
        '''

        return {}


    def _config_arg_parsing(self) -> None:
        '''
        Configure subclass-specific argument parsing.

        Defaults to a noop.
        '''

        pass


    def _show_header(self) -> None:
        '''
        Display a human-readable synopsis of this application, typically by
        logging the basename and current version of this application and various
        metadata assisting debugging of end user issues.

        Defaults to a noop.
        '''

        pass
