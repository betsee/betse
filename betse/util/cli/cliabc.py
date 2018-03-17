#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Top-level abstract base class of all command line interface (CLI) subclasses.
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
from betse import metadata as betse_metadata
from betse import pathtree
from betse.lib import libs
from betse.util.cli.cliarg import SemicolonAwareHelpFormatter
from betse.util.io.log import logs, logconfig
from betse.util.io.log.logenum import LogLevel
from betse.util.path.command import cmds
from betse.util.path.command.cmdexit import SUCCESS, FAILURE_DEFAULT
from betse.util.py.pyprofile import profile_callable, ProfileType
from betse.util.type import types
from betse.util.type.decorator.deccls import abstractproperty
from betse.util.type.text import strs
from betse.util.type.types import (
    type_check,
    ArgParserType,
    MappingType,
    ModuleType,
    SequenceTypes,
    SequenceOrNoneTypes,
)

# ....................{ SUPERCLASS                         }....................
class CLIABC(object, metaclass=ABCMeta):
    '''
    Top-level abstract base class of all command line interface (CLI)
    subclasses, suitable for use by both CLI and GUI front-ends for BETSE.

    This superclass provides _no_ explicit support for subcommands. Only
    concrete subclasses _not_ implementing subcommands should directly subclass
    from this superclass. Concrete subclasses implementing subcommands should
    instead subclass the :class`CLISubcommandableABC` superclass.

    Attributes
    ----------
    _arg_list : list
        List of all passed command-line arguments as unparsed raw strings.
    _arg_parser_top : ArgParserType
        Top-level argument parser parsing all top-level options and subcommands
        passed to this application's external CLI command.
    _args : argparse.Namespace
        :mod:`argparse`-specific container of all passed command-line arguments.
        See "Attributes (_args)" below for further details.
    _exit_status : IntOrNoneTypes
        Exit status with which to exit this application as a byte in the range
        ``[0, 255]``. Defaults to
        :attr:`betse.util.path.command.cmdexit.SUCCESS`. Subclass
        implementations of the :meth:`_do` method may explicitly override this
        default on failure as a CLI-oriented alternative to exception handling.
        Note that, although exit status is typically returned directly by
        callables, doing so here is infeasible due to the :meth:`_do` method API
        already returning profiled objects.
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

        # Default to report success as this application's exit status.
        self._exit_status = SUCCESS

        # For safety, nullify all remaining attributes.
        self._arg_list = None
        self._arg_parser_top = None
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

            # (Re-)initialize all mandatory runtime dependencies of this
            # application *AFTER* parsing and handling all logging-specific CLI
            # options and thus finalizing the logging configuration for the
            # active Python process.
            self._init_app_libs()

            # Run the command-line interface (CLI) defined by this subclass,
            # profiled by the type specified by the "--profile-type" option.
            profile_callable(
                call=self._do,
                profile_type=self._profile_type,
                profile_filename=self._profile_filename,
            )
            # libs.die_unless_runtime_optional('networkx')
            # raise ValueError('Test exception handling.')
        except Exception as exception:
            # Handle this exception.
            self._handle_exception(exception)

            # If this application's exit status is still the default and hence
            # has *NOT* been explicitly overriden by the subclass, replace the
            # default status with failure. If this exception provides a
            # system-specific exit status, this status is used; else, the
            # default failure status (i.e., 1) is used.
            #
            # The Windows-specific "winerror" attribute provided by
            # "WindowsError"-based exceptions is ignored. While more
            # fine-grained than the "errno" attribute, "winerror" values are
            # *ONLY* intended to be used internally rather than reported as an
            # exit status to parent processes.
            if self._exit_status == SUCCESS:
                self._exit_status = getattr(exception, 'errno', FAILURE_DEFAULT)

        # Report this application's exit status to the parent process.
        return self._exit_status

    # ..................{ EXPANDERS                          }..................
    @type_check
    def expand_help(self, text: str, **kwargs) -> str:
        '''
        Interpolate the passed keyword arguments into the passed help string
        template, stripping all prefixing and suffixing whitespace from this
        template.

        For convenience, the following default keyword arguments are
        unconditionally interpolated into this template:

        * ``{script_basename}``, expanding to the basename of the Python wrapper
          script running the current application (e.g., ``betse``).
        * ``{program_name}``, expanding to the human-readable name of this
          application (e.g., ``BETSE``).
        '''

        return strs.remove_whitespace_presuffix(text.format(
            program_name=self._module_metadata.NAME,
            script_basename=cmds.get_current_basename(),
            **kwargs
        ))

    # ..................{ PROPERTIES                         }..................
    @property
    def arg_parser_kwargs(self) -> MappingType:
        '''
        Dictionary of all keyword arguments to be passed to the
        :meth:`ArgParserType.__init__` and
        :meth:`ArgParserType.add_parser` methods of both the top-level argument
        parser *and* all argument subparsers for this CLI.

        Since the :mod:`argparse` API provides multiple methods rather than a
        single method for creating argument parsers, this versatile dictionary
        is preferred over a monolithic factory-based approach (e.g., a
        hypothetical ``_make_arg_parser_top()`` method).
        '''

        return {
            # Wrap non-indented lines in help and description text as paragraphs
            # while preserving indented lines in such text as is.
            'formatter_class': SemicolonAwareHelpFormatter,
        }


    @property
    def arg_parser_top_kwargs(self) -> MappingType:
        '''
        Dictionary of all keyword arguments to be passed to the
        :meth:`ArgP"arserType.__init__` method of the top-level argument parser
        for this CLI.

        See Also
        ----------
        :meth:`_init_arg_parser_top`
            Method initializing this argument parser with this dictionary *and*
            keyword arguments common to all argument parsers.
        '''

        # Dictionary of all keyword arguments to be returned.
        arg_parser_top_kwargs = {
            # Human-readable multi-sentence application description. Since this
            # description is general-purpose rather than CLI-specific, format
            # substrings are *NOT* safely interpolatable into this string.
            'description': self._module_metadata.DESCRIPTION,

            # Human-readable multi-sentence application help suffix.
            'epilog': self.expand_help(self._help_epilog),

            # Basename of the Python wrapper script producing this process.
            'prog': cmds.get_current_basename(),
        }

        # Merge these arguments with all default arguments.
        arg_parser_top_kwargs.update(self.arg_parser_kwargs)

        # Return this dictionary.
        return arg_parser_top_kwargs

    # ..................{ PROPERTIES ~ matplotlib            }..................
    @property
    def _is_option_matplotlib_backend(self) -> bool:
        '''
        ``True`` only if this CLI exposes the ``--matplotlib-backend`` option,
        permitting users to externally specify an arbitrary matplotlib backend at
        the command line.

        Design
        ----------
        Defaults to ``True``. Subclasses overriding this default to ``False``
        should explicitly specify the desired matplotlib backend to use by
        overriding the :meth:`_init_app_libs` method to manually pass the name
        of that backend to the :func:`libs.init`` function
        (e.g., ``libs.init(matplotlib_backend_name='Qt5Agg')``).
        '''

        return True

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
        self._args = self._arg_parser_top.parse_args(self._arg_list)

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

        # Top-level argument parser parsing all top-level options and
        # subcommands passed to the CLI command.
        self._arg_parser_top = ArgParserType(**self.arg_parser_top_kwargs)

        #FIXME: Shift into a new top-level add_arg_parser_options() function of
        #the "cliutil" submodule, passing this function the result of
        #self._make_options_top().

        # For each top-level option, add an argument parsing this option to this
        # argument subparser.
        for option in self._make_options_top():
            option.add(self._arg_parser_top)

    # ..................{ ARGS ~ options                     }..................
    def _make_options_top(self) -> SequenceTypes:
        '''
        Sequence of all :class:`CLIOptionABC` instances defining the top-level
        CLI options accepted by this application.

        For each such option, a corresponding argument is added to the top-level
        argument parser (i.e., :attr:`_arg_parser_top`).

        **Order is significant,** defining the order that the ``--help`` option
        synopsizes these options in. Options omitted here are *not* parsed by
        argument parsers and are thus effectively ignored.

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

        Returns
        ----------
        SequenceTypes
            Sequence of all such :class:`CLIOptionABC` instances.
        '''

        # Avoid circular import dependencies.
        from betse.util.cli.cliopt import (
            CLIOptionArgEnum,
            CLIOptionArgStr,
            CLIOptionBoolTrue,
            CLIOptionVersion,
        )

        # Singleton logging configuration for the current Python process.
        log_config = logconfig.get()

        # Human-readable version specifier suitable for printing to end users.
        version_output = '{} {}'.format(
            cmds.get_current_basename(), self._module_metadata.VERSION)

        # List of all default top-level options to be returned.
        options_top = [
            CLIOptionBoolTrue(
                short_name='-v',
                long_name='--verbose',
                synopsis='print and log all messages verbosely',
            ),

            CLIOptionVersion(
                short_name='-V',
                long_name='--version',
                synopsis='print program version and exit',
                version=version_output,
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
                synopsis=(
                    'type of profiling to perform (defaults to "{default}"):'
                    '\n;* "none", disabling profiling'
                    '\n;* "call", profiling callables (functions, methods)'
                    '\n;* "size", profiling object sizes (requires "pympler")'
                    # ;* "line", profiling code lines (requires "pprofile")
                ),
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
        ]

        # If conditionally exposing the "--matplotlib-backend" option, do so.
        if self._is_option_matplotlib_backend:
            options_top.append(CLIOptionArgStr(
                long_name='--matplotlib-backend',
                synopsis=(
                    'name of matplotlib backend to use '
                    '(see: "betse info")'
                ),
                var_name='matplotlib_backend_name',
                default_value=None,
            ))

        # Return this list.
        return options_top


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

        # Log (and hence display, by default) a human-readable synopsis of this
        # application. Since logging depends on parsing logging operations, this
        # logging is intentionally deferred from the earliest time at which
        # logging could technically be performed (namely, the body of the
        # betse.ignition.ignite() function) to here. The disadvantage of this
        # otherwise sane approach, of course, is that this logging is deferred.
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
        (Re-)initialize this application *before* performing subsequent logic
        assuming this application to have already been initialized.

        Design
        ----------
        Defaults to (re-)initializing all low-level application logic.
        Subclasses may override this method to perform additional
        initialization, in which case this superclass method should still be
        called to properly initialize this application.
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
        self._module_ignition.reinit()


    def _init_app_libs(self) -> None:
        '''
        (Re-)initialize all mandatory runtime dependencies of this application
        *after* parsing and handling all logging-specific CLI options and thus
        finalizing the logging configuration for the active Python process.

        This initialization integrates the custom logging and debugging schemes
        implemented by these dependencies with those implemented by this
        application.

        Design
        ----------
        Defaults to (re-)initializing all mandatory runtime dependencies of
        BETSE. Subclasses overriding this method to perform additional
        initialization must manually call the :meth:`libs.reinit` method to
        initialize these dependencies.
        '''

        # Name of the matplotlib backend explicitly requested by the user if any
        # *OR* "None" otherwise.
        matplotlib_backend_name = None

        # If this CLI exposes the "--matplotlib-backend" option to users, this
        # name is the value of this option.
        if self._is_option_matplotlib_backend:
            matplotlib_backend_name = self._args.matplotlib_backend_name

        # (Re-)initialize all mandatory runtime dependencies.
        libs.reinit(matplotlib_backend_name=matplotlib_backend_name)


    def _show_header(self) -> None:
        '''
        Display a human-readable synopsis of this application, typically by
        logging the basename and current version of this application and various
        metadata assisting debugging of end user issues.
        '''

        # Log this in a manner suitable for downstream applications requiring
        # BETSE as an upstream dependency (e.g., BETSEE).
        logs.log_info(
            'Welcome to <<'
            '{script_name} {script_version} | '
            '{betse_name} {betse_version} | '
            '{betse_codename}'
            '>>.'.format(
                script_name=self._module_metadata.NAME,
                script_version=self._module_metadata.VERSION,
                betse_name=betse_metadata.NAME,
                betse_version=betse_metadata.VERSION,
                betse_codename=betse_metadata.CODENAME,
            ))

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

        On failure, the subclass implementation of this method should either:

        * Raise an exception, in which case this abstract base class implicitly
          logs this exception and report failure as this application's exit status.
        * Explicitly set the :attr:`_exit_status` instance variable to a
          non-zero integer in the range ``[1, 255]`` (e.g.,
          :attr:`betse.util.path.command.cmdexit.FAILURE_DEFAULT`). Note that,
          although exit status is typically returned directly by callables,
          doing so here is infeasible due to this method already returning
          profiled objects.
        '''

        pass

    # ..................{ SUBCLASS ~ mandatory : property    }..................
    # The following properties *MUST* be implemented by subclasses.

    @abstractproperty
    def _help_epilog(self) -> str:
        '''
        Help string template expanded as the **program epilog** (i.e.,
        human-readable string printed after *all* other text in top-level
        application help output).
        '''

        pass


    @abstractproperty
    def _module_ignition(self) -> ModuleType:
        '''
        Imported :mod:`betse.ignition` submodule specific to this application
        (e.g., :mod:`betsee.guiignition` for the BETSEE GUI).

        This property exists principally to support alternate BETSE frontends.
        '''

        pass


    @abstractproperty
    def _module_metadata(self) -> ModuleType:
        '''
        Imported :mod:`betse.metadata` submodule specific to this application
        (e.g., :mod:`betsee.guimetadata` for the BETSEE GUI).

        This property exists principally to support alternate BETSE frontends.
        '''

        pass

    # ..................{ SUBCLASS ~ optional                }..................
    # The following methods may but need *NOT* be implemented by subclasses.

    def _config_arg_parsing(self) -> None:
        '''
        Configure subclass-specific argument parsing.

        Defaults to a noop.
        '''

        pass
