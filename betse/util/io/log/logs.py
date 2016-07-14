#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry
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
# WARNING: To avoid circular import dependencies, import only modules *NOT*
# importing this module at the top-level. Currently, the following modules
# import this module at the top-level and hence *CANNOT* be imported here:
# "betse.util.os.processes".
#
# Since all other modules should *ALWAYS* be able to safely import this module
# at any level, such circularities are best avoided here rather than elsewhere.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

import logging, sys, traceback
from betse.util.type import types
from betse.util.type.types import type_check
from io import StringIO

# ....................{ GETTERS                            }....................
def get(logger_name: str = None) -> logging.Logger:
    '''
    Get the logger with the passed `.`-delimited name, defaulting to the
    basename of the current process (e.g., `betse`) implying the **global
    logger** (i.e., the default application-wide logger).

    This function expects the `LogConfig` class to have been previously
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
def log_debug(message: str, *args, **kwargs) -> None:
    '''
    Log the passed debug message with the root logger, formatted with the passed
    `%`-style positional and keyword arguments.

    This function expects the `LogConfig` class to have been previously
    instantiated, which globally configures logging.
    '''

    logging.debug(message, *args, **kwargs)


@type_check
def log_info(message: str, *args, **kwargs) -> None:
    '''
    Log the passed informational message with the root logger, formatted with
    the passed `%`-style positional and keyword arguments.

    This function expects the `LogConfig` class to have been previously
    instantiated, which globally configures logging.
    '''

    logging.info(message, *args, **kwargs)


@type_check
def log_warning(message: str, *args, **kwargs) -> None:
    '''
    Log the passed warning message with the root logger, formatted with the
    passed `%`-style positional and keyword arguments.

    This function expects the `LogConfig` class to have been previously
    instantiated, which globally configures logging.
    '''

    logging.warning(message, *args, **kwargs)


@type_check
def log_error(message: str, *args, **kwargs) -> None:
    '''
    Log the passed error message with the root logger, formatted with the
    passed `%`-style positional and keyword arguments.

    This function expects the `LogConfig` class to have been previously
    instantiated, which globally configures logging.
    '''

    logging.error(message, *args, **kwargs)


@type_check
def log_exception(exception: Exception) -> None:
    '''
    Log the passed exception with the root logger.
    '''

    try:
        # Avoid circular import dependencies.
        from betse.util.io import stderrs
        from betse.util.io.log import logconfig
        from betse.util.py import identifiers
        from betse.util.type import regexes, strs

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
        _, _, exc_traceback = sys.exc_info()

        #FIXME: The private traceback._iter_chain() function has been removed by
        #Python 3.5. The canonical solution appears to be... wait for it,
        #
        #* Copying the latest Python 3.4 definition of traceback._iter_chain()
        #  into a new _iter_chain() function.
        #* Call our _iter_chain() copy *ONLY* if traceback._iter_chain() is
        #  unavailable.

        # List of 2-tuples "(exception, traceback)" for all parent
        # exceptions of this exception *AND* this exception (in order), where:
        #
        # * "exception" is each exception.
        # * "traceback" is each exception's traceback.
        #
        # Sadly, this list is only gettable via a private module function.
        exc_parents = traceback._iter_chain(
            exception, exc_traceback)

        # String buffer containing a human-readable synopsis of each
        # exception in this chain, unconditionally output to stderr.
        exc_iota_buffer = StringIO()

        # String buffer containing a non-human-readable traceback of each
        # exception in this chain, conditionally logged to the logfile.
        exc_full_buffer = StringIO()

        # Human-readable header prefixing each such buffer.
        buffer_header = 'Exiting prematurely due to fatal error:\n\n'

        # Initialize these buffers to this header.
        exc_iota_buffer.write(buffer_header)
        exc_full_buffer.write(buffer_header)

        # Append each parent exception and that exception's traceback.
        for exc_parent, exc_parent_traceback in (
            exc_parents):
            # If this exception is a string, append this string to the
            # synopsis buffer as is and continue to the next parent. This is
            # an edge case that should *NEVER* happen... but could.
            if types.is_str(exc_parent):
                exc_iota_buffer.write(exc_parent + '\n')
                continue

            # List of traceback lines, excluding this exception's message.
            exc_traceback_lines = traceback.format_exception(
                type(exc_parent),
                exc_parent,
                exc_parent_traceback)
            exc_traceback_lines.pop()

            # List of exception message lines, excluding traceback and hence
            # consisting only of this exception type and original message.
            exc_message_lines = traceback.format_exception_only(
                type(exc_parent), exc_parent)

            # Append this message as is to the traceback buffer *BEFORE*
            # appending a truncation of this message to the message buffer.
            exc_full_buffer.write(
                strs.join(exc_message_lines))
            #print('exception string: '+ exc_message_lines[-1])

            # Split the last line of this message into a non-human-readable
            # exception class and ideally human-readable exception message.
            # If this exception is not None *AND* is convertable without
            # raising exceptions into a string, both format_exception_only()
            # and _format_final_exc_line() guarantee this line to be
            # formatted as follows:
            #     "${exc_class}: ${exc_message}"
            assert types.is_sequence_nonstr_nonempty(
                exc_message_lines), (
                    types.assert_not_sequence_nonstr_nonempty(
                        exc_message_lines, 'Exception message lines'))
            exc_message_match_groups = (
                regexes.get_match_groups_numbered(
                    exc_message_lines[-1],
                    r'^({})(?:\s*|:\s+(.+))$'.format(
                        identifiers.PYTHON_IDENTIFIER_QUALIFIED_REGEX_RAW)))

            # This message is guaranteed to be prefixed by a class name.
            exc_class_name = exc_message_match_groups[0]

            # This message is *NOT* guaranteed to be prefixed by a non-empty
            # message (e.g., assert statements passed no message).
            #
            # If a non-empty message matched, use that.
            exc_message = None
            if exc_message_match_groups[1] is not None:
                exc_message = exc_message_match_groups[1]
            # Else if a debug assertion failed with no explicit message, use
            # the exception context directly detailing this assertion.
            elif exc_class_name == 'AssertionError':
                exc_message = 'Debug assertion failed: {}'.format(
                    # A traceback line typically contains an internal
                    # newline. The substring preceding this newline details
                    # the file and function containing the corresponding
                    # call; the substring following this newline is this
                    # call. Hence, ignore the former.
                    regexes.remove_substrings(
                        exc_traceback_lines[-1], r'^.+\n\s*'))
            # Else, convert this exception's class name into a
            # human-readable message (e.g., from "FileNotFoundError" to
            # "File not found error."). Well, try... at least!
            else:
                exc_message = strs.uppercase_first_char(
                    identifiers.convert_camelcase_to_whitespaced_lowercase(
                        exc_class_name))
            assert types.is_str_nonempty(exc_message), (
                types.assert_not_str_nonempty(
                    exc_message, 'Exception message'))

            # If this class is "KeyError", this message is the single-quoted
            # name of a non-existent key in a dictionary whose access raised
            # this exception. Since that is non-human-readable, wrap this
            # key in human-readable description.
            if exc_class_name == 'KeyError':
                exc_message = 'Dictionary key {} not found.'.format(
                    exc_message)

            # Append this message to the synopsis buffer. For readability,
            # this message is wrapped to the default terminal width and
            # each wrapped line prefixed by indentation.
            exc_iota_buffer.write(strs.wrap(
                text=exc_message, line_prefix='    '))

            # If this exception has a traceback, append this traceback to
            # the traceback but *NOT* synopsis buffer.
            if exc_parent_traceback:
                # Append a traceback header.
                exc_full_buffer.write(
                    '\nTraceback (most recent call last):\n')

                # Append this traceback.
                exc_full_buffer.write(strs.join(
                    # List of lines formatted from this list.
                    traceback.format_list(
                        # List of stack trace entries from this traceback.
                        traceback.extract_tb(
                            exc_parent_traceback))))

        # Append a random error haiku to the traceback buffer... *BECAUSE*!
        exc_full_buffer.write(
            '\n{}'.format(stderrs.get_haiku_random()))

        # String contents of the traceback buffer.
        exc_full = exc_full_buffer.getvalue()

        # Singleton logging configuration for the current Python process.
        log_config = logconfig.get()

        # If logging to a file...
        if log_config.is_logging_file:
            # If verbosity is disabled, output this synopsis to stderr;
            # else, tracebacks containing this synopsis are already
            # output to stderr by logging performeb below.
            if not log_config.is_verbose:
                # Append a reference to the file being logged to to the
                # synopsis buffer.
                exc_iota_buffer.write(
                    '\n\nFor details, see "{}".'.format(log_config.filename))

                # String contents of the synopsis buffer.
                exc_iota = exc_iota_buffer.getvalue()

                # Print this synopsis.
                stderrs.output(exc_iota)

            # Log tracebacks to the debug level and hence *NOT* stderr
            # by default, confining these tracebacks to the logfile.
            # This is a Good Thing (TM). Tracebacks provide more detail
            # than desirable by the typical user.
            log_debug(exc_full)
        # Else, standard file handles are being logged to. In this case,
        # log tracebacks to the error level and hence stderr. Do *NOT*
        # output the synopsis already output in these tracebacks.
        else:
            log_error(exc_full)
    # If this handling raises an exception, catch and print this exception
    # via the standard Python library, guaranteed not to raise exceptions.
    except Exception:
        stderrs.output_exception(
            heading='log_exception() recursively raised exception:\n')
