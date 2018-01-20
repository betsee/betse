#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level logging configuration.
'''

# ....................{ IMPORTS                            }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To avoid circular import dependencies, avoid importing from *ANY*
# application-specific modules at the top-level -- excluding those explicitly
# known *NOT* to import from this module. Since all application-specific modules
# must *ALWAYS* be able to safely import from this module at any level, these
# circularities are best avoided here rather than elsewhere.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

import logging, os, sys
from betse import metadata, pathtree
from betse.util.io.log.logenum import LogLevel
from betse.util.io.log.logfilter import (
    LogFilterThirdPartyDebug, LogFilterMoreThanInfo)
from betse.util.io.log.logformat import LogFormatterWrap
from betse.util.io.log.loghandle import LogHandlerFileRotateSafe
from betse.util.type.types import type_check
from logging import Handler, RootLogger, StreamHandler
from os import path

# ....................{ GLOBALS                            }....................
# See below for utility functions accessing this singleton.
_config = None
'''
Singleton logging configuration for the current Python process.

This configuration provides access to root logger handlers. In particular, this
simplifies modification of logging levels at runtime (e.g., in response to
command-line arguments or configuration file settings).
'''

# ....................{ CONFIG                             }....................
#FIXME: Update docstring to reflect the new default configuration.
class LogConfig(object):
    '''
    BETSE-specific logging configuration.

    This configuration defines sensible default handlers for the root logger,
    which callers may customize (e.g., according to user-defined settings) by
    calling the appropriate getters.

    Caveats
    ----------
    Since this class' :meth:`__init__` method may raise exceptions, this class
    should be instantiated at application startup by an explicit call to the
    module-level :func:`init` function _after_ establishing default exception
    handling. Hence, this class is *not* instantiated at the end of this module.

    Default Settings
    ----------
    All loggers will implicitly propagate messages to the root logger configured
    by this class, whose output will be:

    * Formatted in a timestamped manner detailing the point of origin (e.g.,
      "[2016-04-03 22:02:47] betse ERROR (util.py:50): File not found.").
    * Labelled as the current logger's name, defaulting to `root`. Since this
      is _not_ a terribly descriptive name, callers are encouraged to replace
      this by an application-specific name.
    * Printed to standard error if the logging level for this output is either
      ``WARNING``, ``ERROR``, or ``CRITICAL``.
    * Printed to standard output if the logging level for this output is
      ``INFO``. Together with the prior item, this suggests that output with a
      logging level of ``DEBUG`` will *not* be printed by default.
    * Appended to the user-specific file specified by the
      :func:`pathtree.get_log_default_filename` function, whose:
        * Level defaults to :data:`logger.ALL`. Hence, *all* messages will be
          logged by default, including low-level debug messages. (This is
          helpful for debugging client-side errors.)
      * Contents will be automatically rotated on exceeding a sensible filesize
        (e.g., 16Kb).

    If the default log levels are undesirable, consider subsequently calling
    such logger's `set_level()` method. Since a desired log level is typically
    unavailable until *after* parsing CLI arguments and/or configuration file
    settings *AND* since a logger is required before such level becomes
    available, this function assumes a sane interim default.

    Attributes
    ----------
    _filename : str
        Absolute or relative path of the file logged to by the file handler,
        defaulting to :func:`pathtree.get_log_default_filename`.
    _logger_root : Logger
        Root logger.
    _logger_root_handler_file : Handler
        Root logger handler appending to the current logfile.
    _logger_root_handler_stderr : Handler
        Root logger handler printing to standard error.
    _logger_root_handler_stdout : Handler
        Root logger handler printing to standard output.
    '''

    # ..................{ INITIALIZERS                       }..................
    def __init__(self):
        '''
        Initialize this logging configuration as documented by the class
        docstring.
        '''

        # Avoid circular import dependencies.

        # Initialize the superclass.
        super().__init__()

        # Initialize all non-property attributes to sane defaults. To avoid
        # chicken-and-egg issues, properties should *NOT* be set here.
        self._filename = pathtree.get_log_default_filename()
        self._logger_root = None
        self._logger_root_handler_file = None
        self._logger_root_handler_stderr = None
        self._logger_root_handler_stdout = None

        # Initialize the root logger.
        self._init_logger_root()

        # Initialize root logger handlers *AFTER* the root logger, as the former
        # explicitly add themselves to the latter.
        self._init_logger_root_handler_std()
        self._init_logger_root_handler_file()

        # Redirect all warnings through the logging framewark *AFTER*
        # successfully performing the above initialization.
        logging.captureWarnings(True)


    def _init_logger_root(self) -> None:
        '''
        Initialize the root logger.

        For safety, this function removes all previously initialized handlers
        from this logger.
        '''

        # Root logger.
        self._logger_root = logging.getLogger()

        # For uniqueness, change the name of the root logger to that of our
        # top-level package "betse" from its ambiguous default "root".
        self._logger_root.name = metadata.PACKAGE_NAME

        # Instruct this logger to entertain all log requests, ensuring these
        # requests will be delegated to the handlers defined below. By default,
        # this logger ignores all log requests with level less than "WARNING",
        # preventing handlers from receiving these requests.
        self._logger_root.setLevel(LogLevel.ALL)

        # For safety, remove all existing handlers from the root logger. While
        # this should *NEVER* be the case for conventional BETSE runs, this is
        # usually the case for functional BETSE tests *NOT* parallelized by
        # "xdist" and hence running in the same Python process. For safety,
        # iterate over a shallow copy of the list of handlers to be removed
        # rather than the actual list being modified here.
        for root_handler in list(self._logger_root.handlers):
            self._logger_root.removeHandler(root_handler)


    def _init_logger_root_handler_std(self) -> None:
        '''
        Initialize root logger handlers redirecting log messages to the standard
        stdout and stderr file handles.
        '''

        # Avoid circular import dependencies.
        from betse.util.path.command import cmds

        # Initialize the stdout handler to:
        #
        # * Log only informational messages by default.
        # * Unconditionally ignore all warning and error messages, which the
        #   stderr handler already logs.
        #
        # Sadly, the "StreamHandler" constructor does *NOT* accept the customary
        # "level" attribute accepted by its superclass constructor.
        self._logger_root_handler_stdout = StreamHandler(sys.stdout)
        self._logger_root_handler_stdout.setLevel(LogLevel.INFO)
        self._logger_root_handler_stdout.addFilter(LogFilterMoreThanInfo())

        # Initialize the stderr handler to:
        #
        # * Log only warning and error messages by default.
        # * Unconditionally ignore all informational and debug messages, which
        #   the stdout handler already logs.
        self._logger_root_handler_stderr = StreamHandler(sys.stderr)
        self._logger_root_handler_stderr.setLevel(LogLevel.WARNING)

        # Avoid printing third-party debug messages to the terminal.
        self._logger_root_handler_stdout.addFilter(LogFilterThirdPartyDebug())
        self._logger_root_handler_stderr.addFilter(LogFilterThirdPartyDebug())

        #FIXME: Consider colourizing this format string.

        # Format standard output and error in the conventional way. For a list
        # of all available log record attributes, see:
        #
        #     https://docs.python.org/3/library/logging.html#logrecord-attributes
        #
        # Note that "{{" and "}}" substrings in format() strings escape literal
        # "{" and "}" characters, respectively.
        stream_format = '[{}] {{message}}'.format(cmds.get_current_basename())

        # Formatter for this format.
        stream_formatter = LogFormatterWrap(fmt=stream_format, style='{')

        # Assign these formatters to these handlers.
        self._logger_root_handler_stdout.setFormatter(stream_formatter)
        self._logger_root_handler_stderr.setFormatter(stream_formatter)

        # Register these handlers with the root logger.
        self._logger_root.addHandler(self._logger_root_handler_stdout)
        self._logger_root.addHandler(self._logger_root_handler_stderr)


    def _init_logger_root_handler_file(self) -> None:
        '''
        Initialize the root logger handler appending log messages to the
        currently open logfile file handle.

        This method is designed to be called multiple times, permitting the
        filename associated with this handler to be modified at runtime.
        '''

        # Avoid circular import dependencies.
        from betse.util.path.command import cmds
        from betse.util.type.numeric import ints

        # Absolute or relative path of the directory containing this file.
        file_dirname = path.dirname(self._filename)

        # Minimum level of messages to be log to disk, defaulting to "INFO".
        file_level = LogLevel.INFO

        # If this handler has already been created...
        if self._logger_root_handler_file is not None:
            # Preserve the previously set minimum level of messages to log.
            file_level = self._logger_root_handler_file.level

            # If the root logger has also already been created, remove this
            # handler from this root logger.
            if self._logger_root is not None:
                self._logger_root.removeHandler(self._logger_root_handler_file)

        # If the path of the directory containing this file is non-empty, create
        # this directory if needed. Note this path is empty when this filename
        # is a pure basename (e.g., when the "--log-file=my.log" option is
        # passed).
        #
        # For safety, this directory is created with standard low-level Python
        # functionality rather than our custom higher-level
        # dirs.make_parent_unless_dir() function. The latter logs this creation.
        # Due to the root logger having not been fully configured yet, calling
        # that function here would induce subtle errors or exceptions.
        if file_dirname:
            os.makedirs(file_dirname, exist_ok=True)

        # Root logger file handler, preconfigured as documented above.
        self._logger_root_handler_file = LogHandlerFileRotateSafe(
            filename=self._filename,

            # Append rather than overwrite this file.
            mode='a',

            # Defer opening this file in a just-in-time manner (i.e., until the
            # first call to this handler's emit() method is called to write the
            # first log via this handler). Why? Because (in no particular
            # order):
            #
            # * If the end user requests that *NO* messages be logged to disk
            #   (e.g., by passing the "--log-level=none" option), this file
            #   should *NEVER* be opened and hence created. The simplest means
            #   of doing so is simply to indefinitely defer opening this file.
            # * Doing so slightly reduces the likelihood (but *NOT* eliminate
            #   the possibility) of race conditions between multiple BETSE
            #   processes attempting to concurrently rotate the same logfile.
            delay=True,

            # Encode this file's contents as UTF-8.
            encoding='utf-8',

            # Maximum filesize in bytes at which to rotate this file, equivalent
            # to 1 MB.
            maxBytes=ints.MiB,

            # Maximum number of rotated logfiles to maintain.
            backupCount=8,
        )

        # Initialize this handler's level to the previously established level.
        self._logger_root_handler_file.setLevel(file_level)

        # Prevent third-party debug messages from being logged to disk.
        self._logger_root_handler_file.addFilter(LogFilterThirdPartyDebug())

        # Linux-style logfile format.
        #
        # Note that the "processName" attribute appears to *ALWAYS* expand to
        # "MainProcess", which is not terribly descriptive. Hence, the name of
        # the current process is manually embedded in this format.
        file_format = (
            '[{{asctime}}] '
            '{} {{levelname}} '
            '({{module}}.py:{{funcName}}():{{lineno}}) '
            '<PID {{process}}>:\n'
            '    {{message}}'.format(cmds.get_current_basename()))

        # Format this file according to this format.
        file_formatter = LogFormatterWrap(fmt=file_format, style='{')
        self._logger_root_handler_file.setFormatter(file_formatter)

        # Register this handler with the root logger.
        self._logger_root.addHandler(self._logger_root_handler_file)

    # ..................{ PROPERTIES ~ logger                }..................
    # Read-only properties prohibiting write access to external callers.

    @property
    def logger_root(self) -> RootLogger:
        '''
        **Root logger** (i.e., transitive parent of all other loggers).
        '''

        return self._logger_root

    # ..................{ PROPERTIES ~ handler               }..................
    @property
    def handler_file(self) -> Handler:
        '''
        Root logger handler appending to the current logfile if file logging is
        enabled *or* ``None`` otherwise.
        '''

        return self._logger_root_handler_file


    @property
    def handler_stderr(self) -> Handler:
        '''
        Root logger handler printing to standard error.
        '''

        return self._logger_root_handler_stderr


    @property
    def handler_stdout(self) -> Handler:
        '''
        Root logger handler printing to standard output.
        '''

        return self._logger_root_handler_stdout

    # ..................{ PROPERTIES ~ level                 }..................
    @property
    def file_level(self) -> LogLevel:
        '''
        Minimum level of messages to log to the file handler.
        '''

        return self._logger_root_handler_file.level


    @file_level.setter
    @type_check
    def file_level(self, file_level: LogLevel) -> None:
        '''
        Set the minimum level of messages to log to the file handler.
        '''

        self._logger_root_handler_file.setLevel(file_level)

    # ..................{ PROPERTIES ~ level : verbose       }..................
    @property
    def is_verbose(self) -> bool:
        '''
        ``True`` only if *all* messages are to be unconditionally logged to the
        stdout handler (and hence printed to stdout).

        Equivalently, this method returns ``True`` only if the logging level for
        the stdout handler is :attr:`LogLevel.ALL`.

        Note that this logging level is publicly retrievable by accessing the
        :attr:`handler_stdout.level` property.
        '''

        return self._logger_root_handler_stdout.level == LogLevel.ALL


    @is_verbose.setter
    @type_check
    def is_verbose(self, is_verbose: bool) -> None:
        '''
        Set the verbosity of the stdout handler.

        This method sets this handler's logging level to:

        * If the passed boolean is ``True``, :attr:`LogLevel.ALL` .
        * If the passed boolean is ``False``, :attr:`LogLevel.INFO`.
        '''

        # Convert the passed boolean to a logging level for the stdout handler.
        self._logger_root_handler_stdout.setLevel(
            LogLevel.ALL if is_verbose else LogLevel.INFO)
        # print('handler verbosity: {} ({})'.format(self._logger_root_handler_stdout.level, ALL))

    # ..................{ PROPERTIES ~ path                  }..................
    @property
    def filename(self) -> str:
        '''
        Absolute or relative path of the file logged to by the file handler.
        '''

        return self._filename


    @filename.setter
    @type_check
    def filename(self, filename: str) -> None:
        '''
        Set the absolute or relative path of the file logged to by the file
        handler.

        Due to flaws in the upstream :mod:`logging` API, this method necessarily
        destroys and recreates the current file handler.
        '''

        # If the passed filename is the same as the current filename, avoid
        # unnecessarily destroying and recreating the file handler. This is
        # technically a negligible optimization, but every little bit helps.
        if self._filename == filename:
            return

        # Classify this filename *BEFORE* recreating the file handler, which
        # accesses this variable.
        self._filename = filename

        # Destroy and recreate the file handler.
        self._init_logger_root_handler_file()

# ....................{ INITIALIZERS                       }....................
def init() -> None:
    '''
    Enable the default logging configuration for the active Python process.
    '''

    # Instantiate this singleton global with the requisite defaults.
    # print('Reinitializing logging.')
    global _config
    _config = LogConfig()

# ....................{ GETTERS                            }....................
def get() -> LogConfig:
    '''
    Singleton logging configuration for the active Python process.
    '''

    return _config


def get_metadata():
    '''
    Ordered dictionary synopsizing the current logging configuration.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.mapping.mapcls import OrderedArgsDict

    # Return this dictionary.
    return OrderedArgsDict(
        'file', _config.filename,
        'file level', _config.file_level.name.lower(),
        'verbose', str(_config.is_verbose).lower(),
    )
