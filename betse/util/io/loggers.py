#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

#FIXME: The following blog post provides useful instructions on deserializing
#logging settings from a YAML-formatted configuration file. Leveraging such
#instructions, we could permit users to setup logging as they see fit (e.g.,
#configuring the logging file and level). This is probably the preferable means
#of doing so, rather than providing a gamut of CLI options:
#
#    http://victorlin.me/posts/2012/08/26/good-logging-practice-in-python

#FIXME: It would be great to augment LoggerConfig() with functionality
#permitting the log filename to be explicitly set *AFTER* such object's
#construction. To support this, define a new getter-setter pair of such class
#called simply "filename". Such property's setter should:
#
#* Remove the existing file handler from the root logger.
#* Create a new file handler writing to the passed file.
#* Set such handler on the root logger.
#
#Why not simply pass such filename to such class's __init__()? Because the
#desired log filename will only be available sometime after startup (e.g., after
#parsing command-line arguments and/or configuration files), whereas we would
#prefer that all LoggerConfig() objects come preconfigured with logfile output.
#Why? Because this ensures that, in the event of sufficiently early exceptions,
#there will still exist a logfile for clients to send us.

'''
Low-level logging facilities.

Logging Levels
----------
Logging levels are integer-comparable according to the conventional semantics of
the `<` comparator such that levels assigned smaller integers are more
**inclusive** (i.e., strictly log more messages than) levels assigned larger
integers: e.g.,

    # "DEBUG" is less than and hence more inclusive than "INFO".
    >>> logging.DEBUG < logging.INFO
    True

Logging Hierarchy
----------
Loggers are hierarchically structured according to their `.`-delimited names.
Since the name of the root logger is *always* the empty string, such logger is
*always* the parent of all user-defined loggers. Such hierarchy is an implicit
consequence of logger names and hence requires no manual intervention (e.g., the
root logger `` is implicitly the parent of a user-defined logger `A` is
implicitly the parent of a user-defined logger `A.B`).

By default, logger messages are implicitly propagated up the logger hierarchy
(e.g., messages to logger `A.B` are implicitly progagated to logger `A` are
are implicitly progagated to the root logger). This is a good thing, as it
permits child loggers to consist essentially *only* of a name; assuming all
loggers except the root logger to be unconfigured, messages will be logged
*only* by the root logger.
'''

# ....................{ IMPORTS                            }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To avoid circular import dependencies, import only modules *NOT*
# importing this module at the top-level. Currently, the following modules
# import this module at the top-level and hence *CANNOT* be imported here:
# "betse.util.system.processes".
#
# Since all other modules should *ALWAYS* be able to safely import this module
# at any level, such circularities are best avoided here rather than elsewhere.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
from betse.util.path import files, dirs
from betse.util.type import ints
from logging import StreamHandler
from logging.handlers import RotatingFileHandler
import logging, sys

# ....................{ CONSTANTS ~ int                    }....................
# Originally, we attempted to dynamically copy such constants from the "logging"
# to the current module. Such attempts succeeded in exposing such constants to
# other modules importing this module but *NOT* to this module itself. Hence,
# such constants are manually copied.

DEBUG = logging.DEBUG
'''
Logging level suitable for debugging messages.
'''

INFO = logging.INFO
'''
Logging level suitable for informational messages.
'''

WARNING = logging.WARNING
'''
Logging level suitable for warning messages.
'''

ERROR = logging.ERROR
'''
Logging level suitable for error messages.
'''

CRITICAL = logging.CRITICAL
'''
Logging level suitable for critical messages.
'''

NONE = logging.CRITICAL + 1024
'''
Logging level instructing that no messages be logged.

Since the `logging` module defines no constants encapsulating the concept of
"none" (that is, of logging nothing), this is an ad-hoc constant expected to be
larger than the largest constant defined by such module.
'''

ALL = logging.NOTSET
'''
Logging level instructing that all messages be logged.

Since the `logging` module defines no constants encapsulating the concept of
"all" (that is, of logging everything), this is an ad-hoc constant expected to
be smaller than the smallest constant defined by such module.
'''

# ....................{ GETTERS                            }....................
def get(logger_name: str = None) -> logging.Logger:
    '''
    Get the logger with the passed `.`-delimited name, defaulting to the
    basename of the current process (e.g., `betse`).

    Logger Name
    ----------
    By convention, logger names are typically that of the calling module (e.g.,
    `__name__`). However, since the default logger configuration already
    dynamically embeds the name of such module but *not* such logger into log
    records, there currently remains no reason to pass a logger name.
    '''
    # Default the name of such logger to the basename of the current process.
    # (e.g., "betse").
    if not logger_name:
        from betse.util.system import processes
        logger_name = processes.get_current_basename()

    # If such name is the empty string, this function would get the root logger.
    # Since logging under the root logger is inherently unsafe, assert such
    # constraint.
    assert isinstance(logger_name, str),\
        '"{}" not a string.'.format(logger_name)
    assert logger_name, 'Logger name empty.'

    # Get such logger.
    return logging.getLogger(logger_name)

# ....................{ CONFIG                             }....................
class LoggerConfig(object):
    '''
    Default logging configuration.

    Such configuration defines sensible default handlers for the root logger,
    which callers may customize (e.g., according to user-defined settings) by
    calling the appropriate getters.

    Caveats
    ----------
    Since this class' `__init__()` method may raise exceptions, this class
    should be instantiated at application startup *after* establishing default
    exception handling.

    Default Settings
    ----------
    All loggers will implicitly propagate messages to the root logger configured
    by this class, whose output will be:

    * Formatted in a timestamped manner detailing the point of origin (e.g.,
      "[2016-04-03 22:02:47] betse ERROR (util.py:50): File not found.").
    * Labelled as the current logger's name, defaulting to "root". Since this
      is *not* a terribly descriptive name, callers are encouraged to
    * Printed to standard error if the logging level for such output is either
      `WARNING`, `ERROR`, or `CRITICAL`.
    * Printed to standard output if the logging level for such output is
      `INFO`. (Together with the prior statement, this implies output with a
      logging level of `DEBUG` will *NOT* be printed by default.)
    * Appended to the user-specific file given by `DEFAULT_LOG_FILE`, whose:
      * Level defaults to `logger.ALL`. Hence, *all* messages will be logged by
        default, including low-level debug messages. (This is helpful for
        debugging client-side errors.)
      * Contents will be automatically rotated on exceeding a sensible filesize
        (e.g., 16Kb).

    If the default log levels are undesirable, consider subsequently calling
    such logger's `set_level()` method. Since a desired log level is typically
    unavailable until *after* parsing CLI arguments and/or configuration file
    settings *AND* since a logger is required before such level becomes
    available, this function assumes a sane interim default.

    Attributes
    ----------
    _logger_root : Logger
        Root logger.
    _logger_root_handler_file : Handler
        File handler for the root logger appending to such file.
    _logger_root_handler_stderr : Handler
        Stream handler for the root logger printing to standard error.
    _logger_root_handler_stdout : Handler
        Stream handler for the root logger printing to standard output.
    '''
    def __init__(self):
        '''
        Initialize the root logger for application-wide logging.
        '''
        super().__init__()

        # Import modules required below.
        from betse.util.system import processes

        # Root logger.
        logger_root = logging.getLogger()

        # Instruct such logger to entertain all log requests, ensuring such
        # requests will be delegated to the handlers defined below. By default,
        # such logger ignores all log requests with level less than "WARNING",
        # preventing handlers from receiving such requests.
        logger_root.setLevel(ALL)

        # Root logger stdout handler, preconfigured as documented above. Sadly,
        # such handlers' constructors do *NOT* accept the standard "level"
        # attribute accepted by their base class' constructor.
        self._logger_root_handler_stdout = StreamHandler(sys.stdout)
        self._logger_root_handler_stdout.setLevel(INFO)
        self._logger_root_handler_stdout.addFilter(LoggerFilterInfoOrLess())

        # Root logger stdout handler, preconfigured as documented above.
        self._logger_root_handler_stderr = StreamHandler(sys.stderr)
        self._logger_root_handler_stderr.setLevel(WARNING)

        # If the directory containing such logfile does not exist, fail.
        dirs.die_unless_parent_found(files.DEFAULT_LOG_FILE)

        # Root logger file handler, preconfigured as documented above.
        self._logger_root_handler_file = RotatingFileHandler(
            filename = files.DEFAULT_LOG_FILE,

            # Append rather than overwrite such file.
            mode = 'a',

            # Encode such file's contents as UTF-8.
            encoding = 'utf-8',

            # Filesize at which to rotate such file.
            maxBytes = 16 * ints.KB,

            # Maximum number of rotated logfiles to maintain.
            backupCount = 8,
        )
        self._logger_root_handler_file.setLevel(ALL)

        # Basename of the current process (e.g., "betse").
        script_basename = processes.get_current_basename()

        #FIXME: Colourize me please.

        # Format stdout and stderr output in the conventional way. For a list of
        # all available log record attributes, see:
        #
        #     https://docs.python.org/3/library/logging.html#logrecord-attributes
        #
        # Note that the "processName" attribute appears to *ALWAYS* expand to
        # "MainProcess", which is not terribly descriptive. Hence, the name of
        # the current process is manually embedded in such format.
        #
        # Note that "{{" and "}}" substrings in format() strings escape literal
        # "{" and "}" characters, respectively.
        stream_format = '[{}] {{message}}'.format(script_basename)
        self._logger_root_handler_stdout.setFormatter(logging.Formatter(
            stream_format, style='{',))
        self._logger_root_handler_stderr.setFormatter(logging.Formatter(
            stream_format, style='{',))

        # Enforce a Linux-style logfile format.
        file_format = '[{{asctime}}] {} {{levelname}} ({{module}}.py:{{funcName}}():{{lineno}}):\n    {{message}}'.format(
            script_basename)
        self._logger_root_handler_file.setFormatter(logging.Formatter(
            file_format, style='{',))

        # Register such handlers with the root logger.
        logger_root.addHandler(self._logger_root_handler_stdout)
        logger_root.addHandler(self._logger_root_handler_stderr)
        logger_root.addHandler(self._logger_root_handler_file)

    # ..................{ PROPERTIES                         }..................
    @property
    def file(self) -> logging.Handler:
        '''
        Get the root logger handler appending to the current logfile.
        '''
        return self._logger_root_handler_file

    @property
    def stderr(self) -> logging.Handler:
        '''
        Get the root logger handler printing to standard error.
        '''
        return self._logger_root_handler_stderr

    @property
    def stdout(self) -> logging.Handler:
        '''
        Get the root logger handler printing to standard output.
        '''
        return self._logger_root_handler_stdout

    # ..................{ GETTERS                            }..................
    def get_logger(self, *args, **kwargs) -> logging.Logger:
        '''
        Get the logger with the passed `.`-delimited name, defaulting to the
        basename of the current process (e.g., `betse`).

        This is a convenience method wrapping this module's `get()` function.
        See such function for further details.
        '''
        return get(*args, **kwargs)

# ....................{ CONFIG                             }....................
class LoggerFilterInfoOrLess(logging.Filter):
    '''
    Filter ignoring log records with logging level strictly greater than `INFO`.

    This filter retains only log records with logging level of either `INFO` or
    `DEBUG`.
    '''
    def filter(self, log_record: logging.LogRecord) -> str:
        '''
        True if the passed log record has a logging level of `INFO` or less.
        '''
        assert isinstance(log_record, logging.LogRecord),\
            '"{}" not a log record.'.format(log_record)
        return log_record.levelno <= logging.INFO

# --------------------( WASTELANDS                         )--------------------
        # Prevent the root logger from ignoring *ANY* log requests. (By default,
        # such logger ignores
        # stream_format = '[{processName}] {message}'
        # file_format =\
        #     '[{asctime}] {processName} {levelname} ({module}.py:{funcName}():{lineno}):\n{message}'

            # '[{asctime}] {processName} {levelname:8s} ({module}.py:{funcName}():{lineno}):\n{message}'
# def _copy_logging_levels():
#     '''
#     Copy all logging levels from the `logging` module complete with docstrings
#     into this module on the first importation of this module. Following such
#     importation, all such levels will be publicly accessible as
#     `logger.{logging_level_name}` (e.g., `logger.DEBUG`).
#     '''
#     import sys
#     # thismodule = sys.modules[__name__]
#
#     module_constant = globals()
#     for logging_level_name in ('DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL',):
#         logging_level = getattr(logging, logging_level_name)
#         # setattr(thismodule, logging_level_name, logging_level)
#         globals()[logging_level_name] = logging_level
#
#     global INFO
#     print('INFO: '+INFO)
#
# # Make it so.
# _copy_logging_levels
#
# # For mild safety, delete such function after use. Python's lack of support for
# # multiline anonymous functions becomes unctuous, with time.
# del(_copy_logging_levels)

# class LoggerFilterInfoOnly(logging.Filter):
#     '''
#     Filter ignoring log records with logging level *not* `INFO`.
#
#     This filter retains only log records with a logging level of `INFO`.
#     '''
#     def filter(self, log_record: logging.LogRecord) -> str:
#         '''
#         True if the passed log record has a logging level of `INFO`.
#         '''
#         assert isinstance(log_record, logging.LogRecord),\
#             '"{}" not a string.'.format(log_record)
#         return log_record.levelno == logging.INFO

    # Since logging under the root logger is inherently unsafe, an exception is
    # raised.
    # if logger_name == '':
    #     raise BetseExceptionLog('Logger name empty.')
# from betse.util.exceptions import BetseExceptionLog
        # Note that "{{" and "}}" substrings in format() strings escape literal
        # "{" and "}" characters, respectively.
        # file_format =\
        #     '[{{asctime}}] {} {{levelname:8s}} ({{module}}.py:{{lineno}}): {{message}}'.format(
        #         script_basename)
        # stream_format = '[{}] {{message}}'.format(script_basename)
        # Basename of the current process (e.g., "betse").
        # script_basename = process.get_current_basename()

    # Ideally, this function would automatically get such name by inspecting the
    # current call stack for such name. Unfortunately, a recent stackoverflow
    # comment suggests such inspection fails under PyInstaller: "Also note that
    # this logic fails to work properly when you compile your python code into an
    # exe using pyinstaller." See also:
    #
    #     https://stackoverflow.com/questions/1095543/get-name-of-calling-functions-module-in-python
    # This function should *always* be called on a module-specific basis before
    # attempting to log. In particular, the root logger should *not* be logged to.

            # '[{asctime}] {name} {levelname:8s} ({module}.py:{lineno}): {message}',
        # # Basename of the currently executed external script (e.g., `betse`).
        # script_basename =

    # def __init__(self, script_basename: str):
    #     '''
    #     Initialize the root logger for application-wide logging.
    #
    #     Parameters
    #     ----------
    #     script_basename
    #         Basename of the currently executed external script (e.g., `betse`).
    #     '''
        # assert isinstance(script_basename, str),\
        #     '"{}" not a string.'.format(script_basename)
    # Filter ignoring log records whose logging level strictly greater than
    # `INFO`.
    #
    # Equivalently, this filter retains only log records with a logging level of
    # either `INFO` or `DEBUG`.
    # '''
    # def filter(self, log_record: logging.LogRecord) -> str:
    #     '''
    #     True if the passed log record has a logging level of `INFO` or less
    #     (i.e., is either `INFO` or `DEBUG`).
    #     '''
    #     assert isinstance(log_record, logging.LogRecord),\
    #         '"{}" not a string.'.format(log_record)
    #     return log_record.levelno <= logging.INFO

        # self._logger_root_handler_file.setLevel()
        # Root logger wrapper, preconfigured as documented above.
    # Callers may customize such handlers by calling the appropriate getters.
    # requiring user-specified
    # provides convenient access to such handlers

    # Return such wrapper.
    # return logger_root
    # assert isinstance(script_basename, str), '"{}" not a string.'.format(
    #     script_basename)

# def get_root(script_basename: str) -> None:
      # * Labelled as originating from the external script with the passed
      #   basename (e.g., `betse-qt`).
        # ''.join((
        #     '[{asctime}] {name}',
        #     script_basename,
        #     ' {levelname:8s} ({module}.py:{lineno}): {message}',
        #     '{asctime}   {message}',
        # )),

# such that One logging level being less than another implies the
# former is more inclusive than the latter:
# ....................{ CLASSES                            }....................
# ....................{ SINGLETONS                         }....................
#FUXME: Implement me.
# logger = None
# '''
# Singleton logger usable by both the CLI and GUI interfaces.
# '''
