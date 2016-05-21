#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Abstract command line interface (CLI).
'''

# ....................{ IMPORTS                            }....................
import sys
from abc import ABCMeta, abstractmethod
from argparse import ArgumentParser
from betse import ignition, metadata, pathtree
from betse.cli import info
from betse.util.command import commands
from betse.util.command.args import HelpFormatterParagraph
from betse.util.command.exits import SUCCESS, FAILURE_DEFAULT
from betse.util.io.log import logs
from betse.util.io.log.log_config import LogType, log_config
from betse.util.type import types

# ....................{ CLASS                              }....................
class CLIABC(metaclass=ABCMeta):
    '''
    Abstract command line interface (CLI) suitable for use by both CLI and GUI
    front-ends for BETSE.

    Attributes
    ----------
    _arg_list : list
        List of all passed command-line arguments as unparsed raw strings.
    _arg_parser : ArgumentParser
        `argparse`-specific parser of command-line arguments.
    _arg_parser_kwargs : dict
        Dictionary of keyword arguments which which to initialize
        `ArgumentParser` instances, suitable for passing to both the
        `ArgumentParser.__init__()` and `ArgumentParser.add_parser()` methods.
    _args : argparse.Namespace
        `argparse`-specific object of all passed command-line arguments. See
        "Attributes (_args)" below for further details.
    _script_basename : str
        Basename of the current process (e.g., `betse`).

    Attributes (_args)
    ----------
    is_verbose : bool
        `True` only if low-level debugging messages are to be logged. Defaults
        to `False`.
    log_filename : str
        Absolute or relative path of the file to log to if `log_type` is `file`
        _or_ ignored otherwise. Defaults to the absolute path of the default
        user-specific logfile for `betse` on the current platform.
    log_type : str
        Type of logging if any to be performed. For simplicity, this is a
        `LogType` enumeration member as a lowercase string. Defaults to `file`.
    '''

    # ..................{ INITIALIZERS                       }..................
    def __init__(self):
        super().__init__()

        # Since the basename of the current process is *ALWAYS* available,
        # initialize this basename here for simplicity.
        self._script_basename = commands.get_current_basename()

        # If the active Python interpreter is interpreting a block of arbitrary
        # runtime code passed to this interpreter on the command line via
        # Python's "-c" option (e.g., due to being imported and called as a
        # py.test-based test), this basename is "-c". Convert this non-human-
        # readable string into a human-readable string.
        if self._script_basename == '-c':
            self._script_basename = metadata.SCRIPT_NAME_CLI

        # Initialize these keyword arguments.
        self._arg_parser_kwargs = {
            # Wrap non-indented lines in help and description text as paragraphs
            # while preserving indented lines in such text as is.
            'formatter_class': HelpFormatterParagraph,
        }

        # Initialize these fields to "None" to avoid subtle issues elsewhere.
        self._arg_list = None
        self._arg_parser = None
        self._args = None

    # ..................{ PUBLIC                             }..................
    def run(self, arg_list: list = None) -> int:
        '''
        Run the command-line interface (CLI) defined by the current subclass
        with the passed arguments if non-`None` _or_ with the arguments passed
        on the command line (i.e., `sys.argv`) otherwise.

        Parameters
        ----------
        arg_list : list
            List of zero or more arguments to pass to this interface. Defaults
            to `None`, in which case arguments passed on the command line (i.e.,
            `sys.argv`) will be used instead.

        Returns
        ----------
        int
            Exit status of this interface, in the range `[0, 255]`.
        '''

        # Default unpassed arguments to those passed on the command line,
        # ignoring the first element of "sys.argv" (i.e., the filename of the
        # command from which the current Python process was spawned).
        if arg_list is None:
            arg_list = sys.argv[1:]
        assert types.is_sequence_nonstr(arg_list), (
            types.assert_not_sequence_nonstr(arg_list))
        # print('BETSE arg list (in run): {}'.format(arg_list))

        # Classify arguments for subsequent use.
        self._arg_list = arg_list

        try:
            # Initialize the current application *BEFORE* subsequent logic. This
            # initializes logging and validates paths -- among other chores.
            ignition.init()

            # Parse these arguments *AFTER* initializing logging, ensuring
            # logging of exceptions raised by this parsing.
            self._parse_args()

            # Perform subclass-specific logic.
            self._do()

            # Exit with successful exit status from the current process.
            return SUCCESS
        except Exception as exception:
            # Log this exception.
            logs.log_exception(exception)

            # Exit with failure exit status from the current process. If this
            # exception provides a system-specific exit status, use this status;
            # else, use the default failure status (i.e., 1).
            #
            # Ignore the Windows-specific "winerror" attribute provided by
            # "WindowsError"-based exceptions. While more fine-grained than the
            # "errno" attribute, "winerror" values are *ONLY* intended to be
            # used internally rather than returned as exit status.
            return getattr(exception, 'errno', FAILURE_DEFAULT)

    # ..................{ ARGS                               }..................
    def _parse_args(self) -> None:
        '''
        Parse all currently passed command-line arguments.

        In order, this method:

        * Creates and configures an argument parser with sensible defaults.
        * Calls the subclass-specific `_config_arg_parsing()` method,
          defaulting to a noop.
        * Parses all arguments with such parser.
        '''

        # Configure argument parsing.
        self._make_arg_parser()

        # Parse unnamed string arguments into named, typed arguments.
        self._args = self._arg_parser.parse_args(self._arg_list)

        # Parse top-level options globally applicable to *ALL* subcommands.
        self._parse_global_options()


    def _make_arg_parser(self) -> None:
        '''
        Create and classify the top-level argument parser.
        '''

        # Dictionary of keyword arguments initializing the core argument parser.
        arg_parser_kwargs = {
            # Script name.
            'prog': self._script_basename,

            # Script description.
            'description': metadata.DESCRIPTION,
        }

        # Update this dictionary with preinitialized arguments.
        arg_parser_kwargs.update(self._arg_parser_kwargs)

        # Update this dictionary with subclass-specific arguments.
        arg_parser_kwargs.update(self._get_arg_parser_top_kwargs())

        # Core argument parser.
        self._arg_parser = ArgumentParser(**arg_parser_kwargs)

        # Configure top-level options globally applicable to *ALL* subcommands.
        self._config_global_options()

        # Perform subclass-specific argument parsing configuration.
        self._config_arg_parsing()

    # ..................{ ARGS ~ options                     }..................
    def _config_global_options(self) -> None:
        '''
        Configure argument parsing for top-level options globally applicable to
        _all_ subcommands.
        '''

        # Default values for top-level options configured below, deferred until
        # *AFTER* the ignition.init() function setting these defaults has been
        # called above.
        log_filename_default = pathtree.LOG_DEFAULT_FILENAME
        log_type_default = LogType.FILE.name.lower()

        # Program version specifier.
        program_version = '{} {}'.format(
            self._script_basename, metadata.__version__)

        # Configure top-level options globally applicable to *ALL* subcommands.
        self._arg_parser.add_argument(
            '-v', '--verbose',
            dest='is_verbose',
            action='store_true',
            help='print low-level debugging messages',
        )
        self._arg_parser.add_argument(
            '--log-type',
            dest='log_type',
            action='store',
            choices=tuple(log_type.name.lower() for log_type in LogType),
            default=log_type_default,
            help=(
                'type of logging to perform (defaults to "{}"):\n'
                ';* "none", logging to stdout and stderr\n'
                ';* "file", logging to the "--log-file" file'.format(
                    log_type_default)
            ),
        )
        self._arg_parser.add_argument(
            '--log-file',
            dest='log_filename',
            action='store',
            default=log_filename_default,
            help=(
                'file to log to if "--log-type" is "file" '
                '(defaults to "{}")'.format(log_filename_default)
            ),
        )
        self._arg_parser.add_argument(
            '-V', '--version',
            action='version',
            version=program_version,
            help='print program version and exit',
        )


    def _parse_global_options(self) -> None:
        '''
        Parse top-level options globally applicable to _all_ subcommands.
        '''

        # Configure logging according to the passed options. Note that order of
        # assignment is insignificant here.
        log_config.is_verbose = self._args.is_verbose
        log_config.filename = self._args.log_filename

        # Configure logging type *AFTER* all other logging options. Attempting
        # to set a "log_type" of "FILE" before setting a "filename" will raise
        # an exception, as sanity demands. To do so, this logging type is
        # converted from a lowercase string into an uppercase enumeration
        # member. Since the former is guaranteed by the argument parsing
        # configuration above to be valid, validation need *NOT* be performed.
        log_config.log_type = LogType[self._args.log_type.upper()]

        # Log a one-line synopsis of metadata logged by the "info" subcommand.
        info.log_header()

        # Log all string arguments passed to this command.
        logs.log_debug('Passed argument list {}.'.format(self._arg_list))

    # ..................{ UTILITIES                          }..................
    def _format_help_template(self, text: str) -> str:
        '''
        Format the passed help string template.

        Specifically, this method replaces all instances in this template of:

        * `{script_basename}` by the basename of the current script (e.g.,
          `betse`).
        * `{program_name}` by the name of the current program (e.g., `BETSE`).
        '''
        assert types.is_str(text), types.assert_not_str(text)
        return text.format(
            program_name=metadata.NAME,
            script_basename=self._script_basename,
        )

    # ..................{ SUBCLASS ~ mandatory               }..................
    # The following methods *MUST* be implemented by subclasses.

    @abstractmethod
    def _do(self):
        '''
        Perform subclass-specific logic.
        '''
        pass

    # ..................{ SUBCLASS ~ optional                }..................
    # The following methods may but need *NOT* be implemented by subclasses.

    def _get_arg_parser_top_kwargs(self):
        '''
        Get a subclass-specific dictionary of keyword arguments to be passed to
        the top-level argument parser.
        '''
        return {}


    def _config_arg_parsing(self):
        '''
        Configure subclass-specific argument parsing.
        '''
        pass
