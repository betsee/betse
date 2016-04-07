#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Abstract command line interface (CLI).
'''

# ....................{ IMPORTS                            }....................
from abc import ABCMeta, abstractmethod
from argparse import ArgumentParser
from betse import ignition, metadata
from betse.util.io import loggers, stderr
from betse.util.py import identifiers
from betse.util.os import processes
from betse.util.os.args import HelpFormatterParagraph
from betse.util.type import regexes, strs, types
from io import StringIO
import sys, traceback

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
        Dictionary of keyword arguments to initialize `ArgumentParser` objects
        with. This dictionary is suitable for passing to both the
        `ArgumentParser` constructor *and* the `add_parser()` method of
        `ArgumentParser` objects.
    _args : argparse.Namespace
        `argparse`-specific object of all passed command-line arguments.
    _script_basename : str
        Basename of the current process (e.g., `betse`).
    '''
    def __init__(self):
        super().__init__()

        # Since the basename of the current process is *ALWAYS* available,
        # initialize such basename here for simplicity.
        self._script_basename = processes.get_current_basename()

        # Initialize such keyword arguments.
        self._arg_parser_kwargs = {
            # Wrap non-indented lines in help and description text as paragraphs
            # while preserving indented lines in such text as is.
            'formatter_class': HelpFormatterParagraph,
        }

        # Initialize such fields to None to avoid subtle issues elsewhere (e.g.,
        # attempting to access such logger within _print_exception()).
        self._arg_parser = None
        self._args = None

    # ..................{ PUBLIC                             }..................
    def run(self) -> int:
        '''
        Run the command line interface (CLI) defined by the current subclass.

        Returns
        ----------
        int
            Exit status of such interface, guaranteed to be a non-negative
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
        # Program version specifier.
        program_version = '{} {}'.format(
            self._script_basename, metadata.__version__)

        # Dictionary of keyword arguments with which to initialize the top-level
        # argument parser.
        arg_parser_kwargs = {
            # Script name.
            'prog': self._script_basename,

            # Script description.
            'description': metadata.DESCRIPTION,
        }

        # Update such dictionary with preinitialized such arguments.
        arg_parser_kwargs.update(self._arg_parser_kwargs)

        # Update such dictionary with subclass-specific such arguments.
        arg_parser_kwargs.update(self._get_arg_parser_top_kwargs())

        # Make a command-line argument parser.
        self._arg_parser = ArgumentParser(**arg_parser_kwargs)
        self._arg_parser.add_argument(
            '-v', '--verbose',
            dest = 'is_verbose',
            action = 'store_true',
            help = 'print low-level debugging messages',
        )
        self._arg_parser.add_argument(
            '-V', '--version',
            action = 'version',
            version = program_version,
            help='print program version and exit',
        )

        # Perform subclass-specific argument parsing configuration.
        self._configure_arg_parsing()

        # Parse arguments.
        self._args = self._arg_parser.parse_args()

        # If the user requested verbosity, set the log level for the standard
        # output logger handler to the all-inclusive "ALL".
        if self._args.is_verbose:
            loggers.config.handler_stdout.setLevel(loggers.ALL)

    def _format_help_template(self, text: str) -> str:
        '''
        Format the passed help string template.

        Specifically:

        * Replace all instances in such template of:
          * `{script_basename}` by the basename of the current script (e.g.,
            `betse`).
          * `{program_name}` by the name of the current program (e.g., `BETSE`).
        '''
        assert isinstance(text, str), '"{}" not a string.'.format(text)
        return text.format(
            program_name = metadata.NAME,
            script_basename = self._script_basename,
        )

    # ..................{ EXCEPTIONS                         }..................
    def _print_exception(self, exception: Exception) -> None:
        '''
        Print the passed exception to standard error *and* log such exception.
        '''
        assert isinstance(exception, Exception),\
            '"{}" not an exception.'.format(exception)

        try:
            # Log a descriptive header as an error, thus printing such header to
            # standard error as well by default.
            loggers.log_error(
                'Exiting prematurely due to fatal error:\n')

            # While all loggers provide an exception() method for logging
            # exceptions, the output produced by such method is in the same
            # format as that produced by the Python interpreter on uncaught
            # exceptions. In order, this is:
            #
            # * Such exception's non-layman-readable stack trace.
            # * Such exception's layman-readable error message.
            #
            # Since such format is (arguably) unreadable for non-developers,
            # such exception is reformatted for readability. Sadly, this
            # precludes us from calling our logger's exception() method.

            # Traceback object for such exception.
            _, _, exception_traceback = sys.exc_info()

            # List of tuple pairs comprising both the parent exceptions of
            # such exception *AND* such exception. The first and second
            # items of such pairs are those exceptions and those exception's
            # tracebacks respectively. (Sadly, this is only accessible as a
            # private module function.)
            exception_parents = traceback._iter_chain(
                exception, exception_traceback)

            # String buffer to be printed to standard error, formatting such
            # exception and all parents of such exception in the current
            # exception chain in a user-centric [read: terse] manner.
            stderr_buffer = StringIO()

            # String buffer to be logged to the current logfile, formatting
            # such metadata in a developer-centric [read: verbose] manner.
            log_buffer = StringIO()

            # Append each parent exception and such exception's traceback.
            for exception_parent, exception_parent_traceback in\
                exception_parents:
                # If such exception is a string, append such string to such
                # buffers as is and continue to the next parent.
                if isinstance(exception_parent, str):
                    stderr_buffer.write(exception_parent + '\n')
                    log_buffer   .write(exception_parent + '\n')
                    continue

                # List of traceback lines, excluding the exception message.
                exception_traceback_lines = traceback.format_exception(
                    type(exception_parent),
                    exception_parent,
                    exception_parent_traceback)
                exception_traceback_lines.pop()

                # List of exception message lines, excluding traceback and hence
                # consisting only of this exception type and original message.
                exception_message_lines = traceback.format_exception_only(
                    type(exception_parent), exception_parent)

                # Append this message to the log buffer *BEFORE* appending this
                # message to the standard error buffer. (The latter requires
                # truncating this message for human-readability.)
                log_buffer.write(strs.join(exception_message_lines))
                #print('exception string: '+ exception_message_lines[-1])

                # Split the last line of this message into a non-human-readable
                # exception class and ideally human-readable exception message.
                # If this exception is not None *AND* is convertable without
                # raising exceptions into a string, both format_exception_only()
                # and _format_final_exc_line() guarantee this line to be
                # formatted as follows:
                #     "${exception_class}: ${exception_message}"
                assert len(exception_message_lines),\
                    'Exception message lines empty.'
                exception_message_match_groups = \
                    regexes.get_match_groups_numbered(
                        exception_message_lines[-1],
                        r'^({})(?:\s*|:\s+(.+))$'.format(
                            identifiers.PYTHON_IDENTIFIER_QUALIFIED_REGEX_RAW))

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

                # Append this message to the standard error buffer. For
                # readability, wrap this message to the default terminal width
                # and prefix each wrapped line with indentation.
                stderr_buffer.write(
                    strs.wrap(
                        text = exception_message,
                        line_prefix = '    ',))

                # If such exception has a traceback, append this traceback to
                # this log but *NOT* standard error buffer.
                if exception_parent_traceback:
                    # Append a traceback header.
                    log_buffer.write('\nTraceback (most recent call last):\n')

                    # Append such traceback.
                    log_buffer.write(strs.join(
                        # List of lines formatted from this list.
                        traceback.format_list(
                            # List of stack trace entries from this traceback.
                            traceback.extract_tb(
                                exception_parent_traceback))))

            # Append a logfile reference to this standard error message.
            stderr_buffer.write('\n\nFor details, see "{}".'.format(
                loggers.config.filename))

            # Append a random error haiku to this log message.
            log_buffer.write('\n{}\n'.format(stderr.get_haiku_random()))

            # Exception messages.
            stderr_message = stderr_buffer.getvalue()
            log_message = log_buffer.getvalue()

            # True if the user requested verbosity.
            is_verbose = getattr(self._args, 'is_verbose', False)

            # If a logger has been initialized, log this exception as a debug
            # message. Unless the user explicitly passed command-line option
            # "--verbose" to this script, logging with the debug level confines
            # such traceback to the logfile. This is a (largely) good thing;
            # tracebacks convey more details than expected by customary users.
            if loggers.config.is_initted:
                # Log such message.
                loggers.log_debug(log_message)

                # If the user did *NOT* request verbosity, print a terse message
                # to standard error.
                if not is_verbose:
                    stderr.output(stderr_message)
            # Else, print such exception to standard error. Since the log
            # message is more verbose than and hence subsumes the standard error
            # message, only the former is printed.
            else:
                stderr.output(log_message)
        # If such printing raises an exception, catch and print such exception
        # via the standard Python library, guaranteed not to raise exceptions.
        except Exception:
            stderr.output('_print_exception() recursively raised exception:\n')
            traceback.print_exc()

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

    def _configure_arg_parsing(self):
        '''
        Configure subclass-specific argument parsing.
        '''
        pass
