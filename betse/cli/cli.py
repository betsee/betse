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
from betse.util.system import processes
from betse.util.system.args import HelpFormatterParagraph
from betse.util.type import regexes, strs
from io import StringIO
import sys, traceback

# ....................{ CLASS                              }....................
class CLI(metaclass = ABCMeta):
    '''
    Abstract command line interface (CLI) suitable for use by both the front-
    facing CLI and GUI applications for `betse`.

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
            # Initialize the current application *BEFORE* subsequent logic.
            ignition.init()

            # Parse CLI arguments *AFTER* logging, ensuring that exceptions
            # raised by such parsing will be logged.
            self._parse_args()
            # from six.moves import _dummy_thread
            # from six.moves import tkinter

            # Perform subclass-specific logic.
            self._run()

            # Exit with successful exit status from the current process.
            return 0
        except Exception as exception:
            # Print such exception.
            self._print_exception(exception)

            # Exit with failure exit status from the current process. If such
            # exception provides a system-specific exit status, use such status;
            # else, use the default such status.
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

        #FIXME: Enable the "--config" option. More work than we care to invest,
        #at the moment.

        # Add globally applicable arguments.
        # self._arg_parser.add_argument(
        #     '-c', '--config-file',
        #     default = files.DEFAULT_CONFIG_FILE,
        #     dest = 'config_filename',
        #     help = 'config file to read program settings from')
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

                # List of exception message lines, excluding traceback and hence
                # consisting only of such exception type and original message.
                exception_message_lines = traceback.format_exception_only(
                    type(exception_parent), exception_parent)

                # Append such message to the log buffer *BEFORE* appending such
                # message to the standard error buffer. (The latter requires
                # truncating such message for human-readability.)
                log_buffer.write(strs.join(exception_message_lines))

                #FIXME: If such exception type is "KeyError", the remaining
                #exception message consists only of the offending key and hence
                #is non-human-readable. Correct this by capturing the exception
                #type to group 1, testing such type, and responding accordingly.

                # Strip the non-human-readable exception class from the last
                # line of such message. If such exception is not None *AND* is
                # convertable without raising exceptions to a string, both
                # format_exception_only() and _format_final_exc_line() guarantee
                # such line to be formatted as follows:
                #     "${exception_class}: ${exception_message}"
                assert len(exception_message_lines),\
                    'Exception message lines empty.'
                exception_message_lines[-1] = regexes.remove_substrings(
                    exception_message_lines[-1],
                    '^{}:\s+'.format(
                        regexes.PYTHON_IDENTIFIER_QUALIFIED_REGEX_RAW))

                # Append such message to the standard error buffer. For
                # readability, wrap such message to the default terminal width
                # and prefix each wrapped line with indentation.
                stderr_buffer.write(
                    strs.wrap_lines(
                        lines = exception_message_lines,
                        line_prefix = '    ',))

                # If such exception has a traceback, append such traceback to
                # such log but *NOT* standard error buffer.
                if exception_parent_traceback:
                    # Append a traceback header.
                    log_buffer.write('\nTraceback (most recent call last):\n')

                    # Append such traceback.
                    log_buffer.write(strs.join(
                        # List of lines formatted from such list.
                        traceback.format_list(
                            # List of stack trace entries from such traceback.
                            traceback.extract_tb(
                                exception_parent_traceback))))

            # Append a logfile reference to such standard error message.
            stderr_buffer.write('\n\nFor details, see "{}".'.format(
                loggers.config.filename))

            # Append a random error haiku to such log message.
            log_buffer.write('\n' + stderr.get_haiku_random())

            # Exception messages.
            stderr_message = stderr_buffer.getvalue()
            log_message = log_buffer.getvalue()

            # True if the user requested verbosity.
            is_verbose = getattr(self._args, 'is_verbose', False)

            # If a logger has been initialized, log such exception as a debug
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
    def _run(self):
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

# --------------------( WASTELANDS                         )--------------------
# from betse.cli import help
            # is_verbose = getattr(self._args, 'is_verbose', False)
            # If either the user requested verbosity *OR* no loggers have been
            # initialized, print such exception to standard error. Since the log
            # message is more verbose than and hence subsumes the standard error
            # message, only the former is printed.
            # if is_verbose or not loggers.config.is_initted:
            #     stderr.output(log_message)
            # Help text printed *AFTER* all other output when such script is
            # passed no command-line arguments.
            # epilog = self._format_help_template(help.TEMPLATE_EPILOG),

            # Else, print such exception to standard error. Since the log
            # message is more verbose than and hence subsumes the standard error
            # message, print only the former.
            # else:
            #     stderr.output(log_message)

            # raise Exception('Governments, if they endure, always tend increasingly toward aristocratic forms. No government in history has been known to evade this pattern. And as the aristocracy develops, government tends more and more to act exclusively in the interests of the ruling class -- whether that class be hereditary royalty, oligarchs of financial empires, or entrenched bureaucracy.')
            # Initialize such buffers with descriptive headers.
            # stderr_buffer.write(
            #     'Exiting prematurely due to fatal error:\n\n')
            # log_buffer.write(
            #     'Halting prematurely due to uncaught exception:\n\n')

            # Initialize such buffers with descriptive headers.
            # stderr_buffer.write(
            #     'Exiting prematurely due to fatal error:\n\n')
            # log_buffer.write(
            #     'Halting prematurely due to uncaught exception:\n\n')

                # Exception message concatenated from such lines.
                # exception_message = ''.join(exception_message_lines)

            # are likely to feel comfortable viewing.
            # from being printed to either standard output or error.
            # Assuming such logger retains its default
            # configuration, such exception will be propagated up to the root
            # logger and then handled by the stderr handler.
    # ..................{ LOGGING                            }..................
    #         self._configure_logging()
    # def _configure_logging(self) -> None:
    #     '''
    #     Configure the root logger and obtain an application-wide child logger.
    #     '''
    #     # Configure the root logger.
    #     # Configure the root logger.
    #     self._logger_config = LoggerConfig()
    #
    #     # Create an application-wide child logger.
    #     self._logger = self._logger_config.get_logger()

    # _logger_config : LoggerConfig
    #     Logger configuration, providing access to root logger handlers (e.g.,
    #     for modifying logging levels).
    # _logger : Logger
    #     Logger intended to be used globally (i.e., by *all* classes, functions,
    #     and modules) or None if no such logger has been initialized.
# from betse.util.io.loggers import LoggerConfig
        # self._logger_config = None
        # self._logger = None
#FUXME: Can PyInstaller be made to embed setuptools-specific eggs in the
#executable binaries it produces? If not, we'll probably want to avoid even
#calling die_unless_satisfied_all(), as such function is a noop unless such eggs
#are available at runtime.

# from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
            # Print the default values of options in help output.
    # def _configure_arg_parsing(self, arg_parser: ArgumentParser):
    #     '''
    #     Configure subclass-specific argument parsing with the passed top-level
    #     argument parser.
    #     '''
    #     pass
            #FUXME: We probably don't need this. Excise away. Yay!

                # 'Halting prematurely [read: fatally crashing] due to uncaught exception:\n\n')
            # # Else, print such exception via the standard Python library.
            # else:
            #     traceback.print_exc()
#Parse all passed command-line arguments and run the specified command.
                # usage = '{} <command> [<arg>...]'.format(metadata.SCRIPT_NAME_CLI),
        # Logger intended to be used only by CLI modules or None if no such
        # logger has been initialized.
        # Define global command-line arguments (i.e., arguments applicable to
        # all subparsers, added below).
            # Make a child logger specific to this module.
                # script_basename = self._script_basename)
    # @abstractmethod
    # def _script_basename(self) -> str:
    #     '''
    #     Get the basename of the currently executed external script.
    #     '''
    #     pass

                # # Traceback formatted as a list of lines.
                # exception_traceback_lines = traceback.format_list(
                #     traceback.extract_tb(exception_traceback))
                #
                # # Traceback formatted as a list of lines.
                # exception_traceback_lines = traceback.format_list(
                #     traceback.extract_tb(exception_traceback))
                #
                # # Exception formatted as a string concatenating such exception's
                # # message and traceback in that order.
                # # exception_

                    # reversed(exception_parents):
#FUXME: Improve output on uncatched exceptions by wrapping all CLI operations in
#a try-except block that:
#
#* Catches all exceptions.
#* Improves such output.
#* Exits the current process with status 1.

                # self._logger.exception(exception)
            #FUXME: Print such exception.

            # If a logger has been initialized, log such exception *AFTER*
            # printing such exception to standard error. (Logging is
            # substantially more fragile and hence likely to itself raise
            # further exceptions.)
            # if self._logger:
                #FUXME: Erroneous. This will result in duplicate output.
                # self._logger.exception(exception)

        # return exception.errno if hasattr(exception, 'errno') else 1
#FUXME: Define a new module "betse/dependency.py" performing validation of
#external dependencies, both Python and non-Python. Although we believe "yppy"
#implemented such functionality, google about for the optimum Python 3 solution
#to this presumably commonplace problem.

    # _args : list
    #     List of zero or more arguments passed to such interface (e.g., from the
    #     command line).
        # Initialize such arguments to the current argument list, excluding
        # such list's first item. By cross-platform consent, such item is
        # *ALWAYS* the command name for the current process (e.g., "betse") and
        # hence ignorable.
        # self._args = sys.argv[1:]

#List of zero or more external arguments passed from the command line.
# def main(args = None):
#     CLI().run()
#     '''Run betse`'s command line interface (CLI).
#
#     Parameters
#     ----------
#     args : list, optional
#         List of zero or more arguments passed to such interface (e.g., from the
#         command line) or `None` if called as the entry point in an external
#         script installed by `setuptools`.
#     '''
#     # If called from a setuptools-installed script, copy such arguments from the
#     # argument list excluding the first item of such list. By cross-platform
#     # agreement, such item is *ALWAYS* the command name of the current process
#     # (e.g., "betse") and hence ignorable.
#     if args is None:
#         args = sys.argv[1:]

# if __name__ == '__main__':
#     main()

#FUXME; Configure me for CLI usage. Note that I'm no longer convinced that the
#way we launched "yppy" (e.g., "bin/yppy.bash") was ideal. We really want to do
#the "Pythonic" thing here. ruamel.yaml, for example, installs a Python wrapper
#"/usr/lib/yaml" which (in order):
#
#* Finds an appropriate Python interpreter.
#* Replaces the current process with the result of interpreting
#  "/usr/lib/python-exec/python${PYTHON_VERSION}/yaml". Such file appears to be
#  autogenerated by setuptools at installation time.
#FUXME; Hmm: it looks like we want a new file "betse/__main__.py" resembling:
#    from betse.main import main
#    main()
#This then permits betse to be run as follows:
#    # Yes, you either have to be in the parent directory of the directory
#    # containing such "__main__.py" file *OR* you have to fiddle with
#    # ${PYTHONPATH}.
#    >>> cd ~/py/betse
#    >>> python -m betse
#Naturally, this lends itself well to shell scripting. (Yay!)
#FUXME; Wo! Even nicer. setuptools has implicit support for "__main__.py"-style
#entry points. We just need a "setup.py" resembling:
#    setup(
#        # [...]
#        entry_points={
#            'betse': ['betse = betse.main:main'],
#        },
#    )
#What's sweet about this is that we can define additional separate scripts with
#deeper entry points if we need and or want to.
