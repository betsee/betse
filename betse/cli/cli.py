#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Abstract command line interface (CLI).
'''

# ....................{ IMPORTS                            }....................
import sys, traceback
from abc import ABCMeta, abstractmethod
from argparse import ArgumentParser
from betse import ignition, metadata, pathtree
from betse.util.io import loggers, stderr
from betse.util.py import identifiers
from betse.util.os import processes
from betse.util.os.args import HelpFormatterParagraph
from betse.util.type import regexes, strs, types
from enum import Enum
from io import StringIO

# ....................{ ENUMERATIONS                       }....................
LogType = Enum('LogType', ('NONE', 'FILE'))
'''
Enumeration of all possible types of logging performed by `betse`, corresponding
to the global `--log-type` option configured below.

Attributes
----------
none : enum
    Enumeration member redirecting all logging to standard file handles, in
    which case:
    * All `INFO` and `DEBUG` log messages will be printed to stdout.
    * All `ERROR` and `WARNING` log messages will be printed to stderr.
    * All uncaught exceptions will be printed to stderr.
file : enum
    Enumeration member redirecting all logging to the currently configured
    logfile for `betse`.
'''

# ....................{ CLASS                              }....................
class CLI(metaclass = ABCMeta):
    '''
    Abstract command line interface (CLI) suitable for use by both CLI and GUI
    front-ends for `betse`.

    Attributes
    ----------
    _arg_parser : ArgumentParser
        `argparse`-specific parser of command-line arguments.
    _arg_parser_kwargs : dict
        Dictionary of keyword arguments which which to initialize
        `ArgumentParser` instances, suitable for passing to both the
        `ArgumentParser.__init__()` and `ArgumentParser.add_parser()` methods.
    _args : argparse.Namespace
        `argparse`-specific object of all passed command-line arguments. See
        "Attributes (_args)" below for details.
    _is_log_file : bool
        `True` only if logging to a file (i.e., if `_log_type` is
        `LogType.FILE`).
    _is_verbose : bool
        `True` only if low-level debugging messages are to be logged.
    _log_filename : str
        Absolute or relative path of the file to log to if logging to a file
        _or_ ignored otherwise.
    _log_type : LogType
        Type of logging if any to be performed.
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

    def __init__(self):
        super().__init__()

        # Since the basename of the current process is *ALWAYS* available,
        # initialize this basename here for simplicity.
        self._script_basename = processes.get_current_basename()

        # Initialize these keyword arguments.
        self._arg_parser_kwargs = {
            # Wrap non-indented lines in help and description text as paragraphs
            # while preserving indented lines in such text as is.
            'formatter_class': HelpFormatterParagraph,
        }

        # Initialize these fields to "None" to avoid subtle issues elsewhere
        # (e.g., attempting to access this logger within _print_exception()).
        self._arg_parser = None
        self._args = None
        self._is_verbose = None
        self._log_filename = None
        self._log_type = None

    # ..................{ PUBLIC                             }..................
    def run(self) -> int:
        '''
        Command-line interface (CLI) defined by the current subclass.

        Returns
        ----------
        int
            Exit status of this interface, guaranteed to be a non-negative
            integer in `[0, 255]`, where 0 signifies success and all other
            values failure.
        '''
        try:
            # Initialize the current application *BEFORE* subsequent logic. This
            # initializes logging and validates paths -- among other chores.
            ignition.init()

            # Parse CLI arguments *AFTER* initializing logging, ensuring that
            # exceptions raised by such parsing will be logged.
            self._parse_args()

            # Perform subclass-specific logic.
            self._do()

            # Exit with successful exit status from the current process.
            return 0
        except Exception as exception:
            # Print this exception.
            self._print_exception(exception)

            # Exit with failure exit status from the current process. If this
            # exception provides a system-specific exit status, use this status;
            # else, use the default failure status (i.e., 1).
            #
            # Ignore the Windows-specific "winerror" attribute provided by
            # "WindowsError"-based exceptions. While more fine-grained than the
            # "errno" attribute, "winerror" values are *ONLY* intended to be
            # used internally rather than returned as exit status.
            return getattr(exception, 'errno', 1)

    # ..................{ ARGS                               }..................
    def _parse_args(self) -> None:
        '''
        Parse all currently passed command-line arguments.

        In order, this method:

        * Creates and configures an argument parser with sensible defaults.
        * Calls the subclass-specific `_configure_arg_parsing()` method,
          defaulting to a noop.
        * Parses all arguments with such parser.
        '''

        # Configure argument parsing.
        self._make_arg_parser()

        # Parse arguments.
        self._args = self._arg_parser.parse_args()

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
        self._configure_arg_parsing()

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

        #FIXME: Actually use "self._log_filename", please.

        # Classify options requiring no conversion as is.
        self._is_verbose = self._args.is_verbose
        self._log_filename = self._args.log_filename

        #FIXME: Actually use this. As low-hanging fruit, let's support by adding
        #support for this to the _print_exception() method defined below.
        #FIXME: After adding that, we'll probably want to add an additional
        #"log type:" field to "betse info" output. See the "info.py" module.
        #FIXME: After adding that, we'll want to refactor logging to:
        #
        #* Default to "LogType.NONE" and hence to *NOT* log to any files. The
        #  reason why, of course, is that the file to be logged to should be
        #  configurable at runtime. Since chicken-and-the-egg issues rapidly
        #  ensue, the only sane solution is to log only to the terminal on
        #  initial startup.
        #* Add a new method to "betse.util.io.logging", permitting logging to be
        #  reconfigured to log to a file. Shouldn't be terribly arduous.

        # Convert the logging type from a lowercase string into an uppercase
        # enumeration member. Since the former is guaranteed by the
        # configuration above to be valid, validation need *NOT* be performed.
        self._log_type = LogType[self._args.log_type.upper()]

        # True only if logging to a file.
        self._is_log_file = self._log_type is LogType.FILE

        #FIXME: Reconfigure logging to do so, as detailed above.

        # If logging to a file...
        if self._is_log_file:
            pass

        # If the user requested verbosity, set the log level for the standard
        # output logger handler to the all-inclusive "ALL".
        if self._is_verbose:
            loggers.config.handler_stdout.setLevel(loggers.ALL)

    # ..................{ EXCEPTIONS                         }..................
    def _print_exception(self, exception: Exception) -> None:
        '''
        Print the passed exception to standard error *and* log such exception.
        '''
        assert types.is_exception(exception), (
            types.assert_not_exception(exception))

        try:
            # While all loggers provide an exception() method for logging
            # exceptions, the output produced by these methods is in the same
            # format as that produced by the Python interpreter on uncaught
            # exceptions. In order, this is:
            #
            # * This exception's non-human-readable stack trace.
            # * This exception's human-readable error message.
            #
            # Since this format is (arguably) unreadable for non-developers,
            # this exception is reformatted for readability. Sadly, this
            # precludes us from calling our logger's exception() method.

            # Traceback object for this exception.
            _, _, exception_traceback = sys.exc_info()

            # List of 2-tuples "(exception, traceback)" for all parent
            # exceptions of this exception *AND* this exception (in order),
            # where:
            #
            # * "exception" is each exception.
            # * "traceback" is each exception's traceback.
            #
            # Sadly, this list is only gettable via a private module function.
            exception_parents = traceback._iter_chain(
                exception, exception_traceback)

            # String buffer containing a human-readable synopsis of each
            # exception in this chain, unconditionally output to stderr.
            exception_iota_buffer = StringIO()

            # String buffer containing a non-human-readable traceback of each
            # exception in this chain, conditionally logged to the logfile.
            exception_full_buffer = StringIO()

            # Human-readable header prefixing each such buffer.
            buffer_header = 'Exiting prematurely due to fatal error:\n\n'

            # Initialize these buffers to this header.
            exception_iota_buffer.write(buffer_header)
            exception_full_buffer.write(buffer_header)

            # Append each parent exception and that exception's traceback.
            for exception_parent, exception_parent_traceback in (
                exception_parents):
                # If this exception is a string, append this string to the
                # synopsis buffer as is and continue to the next parent. This is
                # an edge case that should *NEVER* happen... but could.
                if types.is_str(exception_parent):
                    exception_iota_buffer.write(exception_parent + '\n')
                    continue

                # List of traceback lines, excluding this exception's message.
                exception_traceback_lines = traceback.format_exception(
                    type(exception_parent),
                    exception_parent,
                    exception_parent_traceback)
                exception_traceback_lines.pop()

                # List of exception message lines, excluding traceback and hence
                # consisting only of this exception type and original message.
                exception_message_lines = traceback.format_exception_only(
                    type(exception_parent), exception_parent)

                # Append this message as is to the traceback buffer *BEFORE*
                # appending a truncation of this message to the message buffer.
                exception_full_buffer.write(
                    strs.join(exception_message_lines))
                #print('exception string: '+ exception_message_lines[-1])

                # Split the last line of this message into a non-human-readable
                # exception class and ideally human-readable exception message.
                # If this exception is not None *AND* is convertable without
                # raising exceptions into a string, both format_exception_only()
                # and _format_final_exc_line() guarantee this line to be
                # formatted as follows:
                #     "${exception_class}: ${exception_message}"
                assert types.is_sequence_nonstr_nonempty(
                    exception_message_lines), (
                        types.assert_not_sequence_nonstr_nonempty(
                            exception_message_lines, 'Exception message lines'))
                exception_message_match_groups = (
                    regexes.get_match_groups_numbered(
                        exception_message_lines[-1],
                        r'^({})(?:\s*|:\s+(.+))$'.format(
                            identifiers.PYTHON_IDENTIFIER_QUALIFIED_REGEX_RAW)))

                # This message is guaranteed to be prefixed by a class name.
                exception_class_name = exception_message_match_groups[0]

                # This message is *NOT* guaranteed to be prefixed by a non-empty
                # message (e.g., assert statements passed no message).
                #
                # If a non-empty message matched, use that.
                exception_message = None
                if exception_message_match_groups[1] is not None:
                    exception_message = exception_message_match_groups[1]
                # Else if a debug assertion failed with no explicit message, use
                # the exception context directly detailing this assertion.
                elif exception_class_name == 'AssertionError':
                    exception_message = 'Debug assertion failed: {}'.format(
                        # A traceback line typically contains an internal
                        # newline. The substring preceding this newline details
                        # the file and function containing the corresponding
                        # call; the substring following this newline is this
                        # call. Hence, ignore the former.
                        regexes.remove_substrings(
                            exception_traceback_lines[-1], r'^.+\n\s*'))
                # Else, convert this exception's class name into a
                # human-readable message (e.g., from "FileNotFoundError" to
                # "File not found error."). Well, try... at least!
                else:
                    exception_message = strs.uppercase_first_char(
                        identifiers.convert_camelcase_to_whitespaced_lowercase(
                            exception_class_name))
                assert types.is_str_nonempty(exception_message), (
                    types.assert_not_str_nonempty(
                        exception_message, 'Exception message'))

                # If this class is "KeyError", this message is the single-quoted
                # name of a non-existent key in a dictionary whose access raised
                # this exception. Since that is non-human-readable, wrap this
                # key in human-readable description.
                if exception_class_name == 'KeyError':
                    exception_message = 'Dictionary key {} not found.'.format(
                        exception_message)

                # Append this message to the synopsis buffer. For readability,
                # this message is wrapped to the default terminal width and
                # each wrapped line prefixed by indentation.
                exception_iota_buffer.write(strs.wrap(
                    text=exception_message, line_prefix='    '))

                # If this exception has a traceback, append this traceback to
                # the traceback but *NOT* synopsis buffer.
                if exception_parent_traceback:
                    # Append a traceback header.
                    exception_full_buffer.write(
                        '\nTraceback (most recent call last):\n')

                    # Append this traceback.
                    exception_full_buffer.write(strs.join(
                        # List of lines formatted from this list.
                        traceback.format_list(
                            # List of stack trace entries from this traceback.
                            traceback.extract_tb(
                                exception_parent_traceback))))

            # Append a random error haiku to the traceback buffer... *BECAUSE*!
            exception_full_buffer.write(
                '\n{}'.format(stderr.get_haiku_random()))

            # If logging to a file, append a reference to this file to the
            # synopsis buffer.
            if self._is_log_file:
                exception_iota_buffer.write(
                    '\n\nFor details, see "{}".'.format(
                        loggers.config.filename))

            # String contents of these buffers.
            exception_iota = exception_iota_buffer.getvalue()
            exception_full = exception_full_buffer.getvalue()

            # If logging has been initialized...
            if loggers.config.is_initted:
                # If logging to a file...
                if self._is_log_file:
                    # If verbosity is disabled, output this synopsis to stderr;
                    # else, tracebacks containing this synopsis are already
                    # output to stderr by logging performeb below.
                    if not self._is_verbose:
                        stderr.output(exception_iota)

                    # Log tracebacks to the debug level and hence *NOT* stderr
                    # by default, confining these tracebacks to the logfile.
                    # This is a Good Thing (TM). Tracebacks provide more detail
                    # than desirable by the typical user.
                    loggers.log_debug(exception_full)
                # Else, standard file handles are being logged to. In this case,
                # log tracebacks to the error level and hence stderr. Do *NOT*
                # output the synopsis already output in these tracebacks.
                else:
                    loggers.log_error(exception_full)
            # Else, print this synopsis and tracebacks directly to stderr.
            else:
                stderr.output(exception_iota)
                stderr.output(exception_full)
        # If this handling raises an exception, catch and print this exception
        # via the standard Python library, guaranteed not to raise exceptions.
        except Exception:
            stderr.output('_print_exception() recursively raised exception:\n')
            traceback.print_exc()

    # ..................{ UTILITIES                          }..................
    def _format_help_template(self, text: str) -> str:
        '''
        Format the passed help string template.

        Specifically, this method replaces all instances in this template of:

        * `{script_basename}` by the basename of the current script (e.g.,
          `betse`).
        * `{program_name}` by the name of the current program (e.g., `BETSE`).
        '''
        assert isinstance(text, str), '"{}" not a string.'.format(text)
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


    #FIXME: Rename to _config_arg_parsing() and likewise for similar subclass
    #methods (e.g., _configure_arg_parsing_plot()).
    def _configure_arg_parsing(self):
        '''
        Configure subclass-specific argument parsing.
        '''
        pass
