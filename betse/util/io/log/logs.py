#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

#FIXME: Error messages should be prefixed by strings uniquely identifying the
#sources of those messages. Specifically:
#
#* Warnings should be prefixed by "<Warning> {source_module_basename}".
#* Errors should be prefixed by "<Error> {source_module_basename}".
#
#See the following stackoveflow question for details on how to implement this:
#    https://stackoverflow.com/questions/14844970/modifying-logging-message-format-based-on-message-logging-level-in-python3

#FIXME: The following blog post provides useful instructions on deserializing
#logging settings from a YAML-formatted configuration file. Leveraging such
#instructions, we could permit users to setup logging as they see fit (e.g.,
#configuring the logging file and level). This is probably the preferable means
#of doing so, rather than providing a gamut of CLI options:
#
#    http://victorlin.me/posts/2012/08/26/good-logging-practice-in-python

'''
Low-level logging facilities.

Logging Hierarchy
----------
Loggers are hierarchically structured according to their `.`-delimited names.
Since the name of the root logger is _always_ the empty string, this logger is
_always_ the parent of all user-defined loggers. This hierarchy is an implicit
consequence of logger names and hence requires no manual intervention (e.g., the
root logger `` is implicitly the parent of a user-defined logger `A` is
implicitly the parent of a user-defined logger `A.B`).

By default, logger messages are implicitly propagated up the logger hierarchy
(e.g., messages to logger `A.B` are implicitly progagated to logger `A` are
are implicitly progagated to the root logger). This is a good thing, effectively
reducing child loggers to symbolic names; assuming all loggers except the root
logger to be unconfigured, messages will be logged _only_ by the root logger.
'''

# ....................{ IMPORTS                            }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To avoid circular import dependencies, avoid importing from *ANY*
# application-specific modules at the top-level -- excluding those explicitly
# known *NOT* to import from this module. Since all application-specific modules
# must *ALWAYS* be able to safely import from this module at any level, these
# circularities are best avoided here rather than elsewhere.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

import logging, sys, traceback
from betse.util.io.log.logenum import LogLevel
from betse.util.type import types
from betse.util.type.types import type_check

# ....................{ GETTERS                            }....................
def get(logger_name: str = None) -> logging.Logger:
    '''
    Get the logger with the passed `.`-delimited name, defaulting to the
    basename of the current process (e.g., `betse`) implying the **global
    logger** (i.e., the default application-wide logger).

    This function expects the :class:`LogConfig` class to have been previously
    instantiated, which globally configures logging.

    Logger Name
    ----------
    By convention, logger names are typically that of the calling module (e.g.,
    `__name__`). However, since the default logger configuration already
    dynamically embeds the name of such module but *not* such logger into log
    records, there currently remains no reason to pass a logger name.
    '''

    # Default the name of this logger to the name of the root logger.
    if logger_name is None:
        logger_name = logging.root.name

    # If this name is the empty string, this function would get the root logger.
    # Since this name being empty typically constitutes an implicit error rather
    # than an attempt to get the root logger, prevent this.
    assert types.is_str_nonempty(logger_name), (
        types.assert_not_str_nonempty(logger_name, 'Logger name'))

    # Return this logger.
    return logging.getLogger(logger_name)

# ....................{ LOGGERS                            }....................
@type_check
def log_levelled(message: str, level: LogLevel, *args, **kwargs) -> None:
    '''
    Log the passed message of the passed logging level (e.g.,
    :attr:`LogLevel.INFO`) with the root logger, formatted with the passed
    ``%``-style positional and keyword arguments.

    This function expects the :class:`LogConfig` class globally configuring
    logging to be instantiated as a singleton.

    Parameters
    ----------
    message : str
        Message to log containing zero or more ``%``-style format specifiers.
    level : LogLevel
        Logging level to log this message with (e.g., :attr:`LogLevel.INFO`).

    ALl remaining parameters are interpolated into the message according to the
    ``%``-style format specifiers embedded in this message.
    '''

    # The Logger.log() method accepts these parameters in the opposite order.
    logging.log(level, message, *args, **kwargs)

# ....................{ LOGGERS ~ level                    }....................
@type_check
def log_debug(message: str, *args, **kwargs) -> None:
    '''
    Log the passed debug message with the root logger, formatted with the passed
    ``%``-style positional and keyword arguments.

    This function expects the :class:`LogConfig` class globally configuring
    logging to be instantiated as a singleton.
    '''

    logging.debug(message, *args, **kwargs)


@type_check
def log_info(message: str, *args, **kwargs) -> None:
    '''
    Log the passed informational message with the root logger, formatted with
    the passed ``%``-style positional and keyword arguments.

    This function expects the :class:`LogConfig` class globally configuring
    logging to be instantiated as a singleton.
    '''

    logging.info(message, *args, **kwargs)


@type_check
def log_warning(message: str, *args, **kwargs) -> None:
    '''
    Log the passed warning message with the root logger, formatted with the
    passed ``%``-style positional and keyword arguments.

    This function expects the :class:`LogConfig` class globally configuring
    logging to be instantiated as a singleton.
    '''

    logging.warning(message, *args, **kwargs)


@type_check
def log_error(message: str, *args, **kwargs) -> None:
    '''
    Log the passed error message with the root logger, formatted with the
    passed ``%``-style positional and keyword arguments.

    This function expects the :class:`LogConfig` class globally configuring
    logging to be instantiated as a singleton.
    '''

    logging.error(message, *args, **kwargs)

# ....................{ LOGGERS ~ exception                }....................
@type_check
def log_exception(exception: Exception) -> None:
    '''
    Log the passed exception with the root logger.
    '''

    # While all loggers provide an exception() method for logging exceptions,
    # the output produced by these methods is in the same format as that
    # produced by the Python interpreter on uncaught exceptions. In order, this
    # is:
    #
    # * This exception's non-human-readable stack trace.
    # * This exception's human-readable error message.
    #
    # Since this format is (arguably) unreadable for non-developers, this
    # exception is reformatted for readability. Sadly, this precludes calling
    # our logger's exception() method.
    #
    # Attempt to...
    try:
        # Avoid circular import dependencies.
        from betse.util.io import exceptions, stderrs
        from betse.util.io.log import logconfig
        from betse.util.io.log.logenum import LogLevel

        # Terse synopsis and verbose traceback for this exception.
        exc_synopsis, exc_traceback = exceptions.get_metadata(exception)

        # Singleton logging configuration for the current Python process.
        log_config = logconfig.get()

        # If the end user requested that nothing be logged to disk, respect this
        # request by logging tracebacks to the error level and hence stderr.
        # (Avoid printing the synopsis already embedded in these tracebacks.)
        if log_config.file_level >= LogLevel.NONE:
            log_error(exc_traceback)
        # Else, the end user requested that at least something be logged to
        # disk. For debuggability, the logging level of the file handler is
        # temporarily decreased to the debug level, guaranteeing that tracebacks
        # are *ALWAYS* at least logged to disk rather than (possibly) discarded.
        # For readability, tracebacks are only logged to stderr if explicitly
        # requested by the end user.
        else:
            # If verbosity is disabled, output this synopsis to stderr;
            # else, tracebacks containing this synopsis are already
            # output to stderr by logging performed below.
            if not log_config.is_verbose:
                # Print this synopsis followed by a human-readable reference to
                # the current logfile.
                stderrs.output(
                    '{}\n\nFor details, see "{}".'.format(
                        exc_synopsis, log_config.filename))

            # Previous minimum level of messages to log to disk.
            log_config_file_level = log_config.file_level

            # Temporarily coerce this to the debug level, ensuring that
            # tracebacks are *ALWAYS* at least logged to disk.
            log_config.file_level = LogLevel.DEBUG

            # Attempt to...
            try:
                # Log tracebacks to the debug level and hence *NOT* stderr by
                # default, isolating tracebacks to disk. This is a Good Thing.
                # Tracebacks supply more detail than desired by typical users.
                log_debug(exc_traceback)
            # Revert to the previous level even if an exception is raised.
            finally:
                log_config.file_level = log_config_file_level
    # If this handling raises an exception, catch and print this exception
    # via the standard Python library, guaranteed not to raise exceptions.
    except Exception:
        # Header preceding the exception to be printed.
        exc_heading = 'log_exception() recursively raised exception:\n'

        # If the stderrs.output_exception() function exists, call that function.
        if 'stderrs' in locals():
            stderrs.output_exception(heading=exc_heading)
        # Else, something has gone horribly wrong. Defer to stock functionality
        # in the standard "traceback" module.
        else:
            print(exc_heading, file=sys.stderr)
            traceback.print_exc(file=sys.stderr)
