#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

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

#FIXME: The following blog post provides useful instructions on deserializing
#logging settings from a YAML-formatted configuration file. Leveraging such
#instructions, we could permit users to setup logging as they see fit (e.g.,
#configuring the logging file and level). This is probably the preferable means
#of doing so, rather than providing a gamut of CLI options:
#
#    http://victorlin.me/posts/2012/08/26/good-logging-practice-in-python

# ....................{ IMPORTS                            }....................
from betse.util.path import dirs
from betse.util.type import ints
import logging, sys

# ....................{ CONSTANTS                          }....................
LOG_FILE = dirs.DOT_DIR
'''
Absolute path of the user-specific file to which `betse` logs messages.
'''

# ....................{ CONSTANTS ~ int                    }....................
def _copy_logging_levels():
    '''
    Copy all logging levels from the `logging` module complete with docstrings
    into this module on the first importation of this module. Following such
    importation, all such levels will be publicly accessible as
    `logger.{logging_level_name}` (e.g., `logger.DEBUG`).
    '''
    module_constant = globals()
    for logging_level_name in ('DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL',):
        logging_level = getattr(logging, logging_level_name)
        module_constant[logging_level_name] = logging_level

# Make it so.
_copy_logging_levels

# For mild safety, delete such function after use. Python's lack of support for
# multiline anonymous functions becomes unctuous, with time.
del(_copy_logging_levels)

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
def get(logger_name: str) -> logging.Logger:
    '''
    Get a new logger with the passed `.`-delimited name.

    This function should *always* be called on a module-specific basis before
    attempting to log. In particular, the root logger should *not* be logged to.

    Logger Name
    ----------
    For simplicity, such name should typically be that of the calling module
    (e.g., `__name__`).

    Ideally, this function would automatically get such name by inspecting the
    current call stack for such name. Unfortunately, a recent stackoverflow
    comment suggests such inspection fails under PyInstaller: "Also note that
    this logic fails to work properly when you compile your python code into an
    exe using pyinstaller." See also:

        https://stackoverflow.com/questions/1095543/get-name-of-calling-functions-module-in-python
    '''
    assert isinstance(logger_name, str),\
        '"{}" not a string.'.format(logger_name)
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
    * Appended to the user-specific file given by `LOG_FILE`, whose:
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
    _logger_root_handler_file : LoggerHandler
        File handler for the root logger appending to such file.
    _logger_root_handler_stderr : LoggerHandler
        Stream handler for the root logger printing to standard error.
    _logger_root_handler_stdout : LoggerHandler
        Stream handler for the root logger printing to standard output.
    '''
    def __init__(self, script_basename: str):
        '''
        Initialize the root logger for application-wide logging.

        Parameters
        ----------
        script_basename
            Basename of the currently executed external script (e.g., `betse`).
        '''
        super().__init__()
        assert isinstance(script_basename, str),\
            '"{}" not a string.'.format(script_basename)

        # If the directory containing such logfile does not exist, fail.
        dirs.die_unless_parent_found(LOG_FILE)

        # Root logger.
        logger_root = logging.getLogger()

        # Root logger stdout handler, preconfigured as documented above.
        self._logger_root_handler_stdout = logging.StreamHandler(
            level = INFO,
            stream = sys.stdout,
        )
        self._logger_root_handler_stdout.addFilter(LoggerFilterInfoOnly())

        # Root logger stdout handler, preconfigured as documented above.
        self._logger_root_handler_stderr = logging.StreamHandler(
            level = WARNING,
            stream = sys.stderr,
        )

        # Root logger file handler, preconfigured as documented above.
        self._logger_root_handler_file = logging.handlers.RotatingFileHandler(
            filename = LOG_FILE,
            level = ALL,

            # Append rather than overwrite such file.
            mode = 'a',

            # Encode such file's contents as UTF-8.
            encoding = 'utf-8',

            # Filesize at which to rotate such file.
            maxBytes = 16 * ints.KB,

            # Maximum number of rotated logfiles to maintain.
            backupCount = 8,
        )

        #FIXME: Colourize me please.

        # Format stdout and stderr output in the conventional way.
        stream_format = ''.join(('[', script_basename, '] {message}',))
        self._logger_root_handler_stdout.setFormatter(logging.Formatter(
            stream_format, style='{',))
        self._logger_root_handler_stderr.setFormatter(logging.Formatter(
            stream_format, style='{',))

        # Enforce a Linux-style logfile format.
        self._logger_root_handler_file.setFormatter(logging.Formatter(
            '[{asctime}] {name} {levelname:8s} ({module}.py:{lineno}): {message}',
            style='{',
        ))

        # Register such handlers with the root logger.
        logger_root.addHandler(self._logger_root_handler_stdout)
        logger_root.addHandler(self._logger_root_handler_stderr)
        logger_root.addHandler(self._logger_root_handler_file)

# ....................{ CONFIG                             }....................
class LoggerFilterInfoOnly(logging.Filter):
    '''
    Filter ignoring log records whose logging level is *not* `INFO`.

    This filter retains only log records with a logging level of `INFO`.
    '''
    def filter(self, log_record: logging.LogRecord) -> str:
        '''
        True if the passed log record has a logging level of `INFO`.
        '''
        assert isinstance(log_record, logging.LogRecord),\
            '"{}" not a string.'.format(log_record)
        return log_record.levelno == logging.INFO

# --------------------( WASTELANDS                         )--------------------
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
