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
from betse.util.type import types
from betse.util.type.types import type_check
from io import StringIO

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
def log_debug(message: str, *args, **kwargs) -> None:
    '''
    Log the passed debug message with the root logger, formatted with the passed
    `%`-style positional and keyword arguments.

    This function expects the :class:`LogConfig` class to have been previously
    instantiated, which globally configures logging.
    '''

    logging.debug(message, *args, **kwargs)


@type_check
def log_info(message: str, *args, **kwargs) -> None:
    '''
    Log the passed informational message with the root logger, formatted with
    the passed `%`-style positional and keyword arguments.

    This function expects the :class:`LogConfig` class to have been previously
    instantiated, which globally configures logging.
    '''

    logging.info(message, *args, **kwargs)


@type_check
def log_warning(message: str, *args, **kwargs) -> None:
    '''
    Log the passed warning message with the root logger, formatted with the
    passed `%`-style positional and keyword arguments.

    This function expects the :class:`LogConfig` class to have been previously
    instantiated, which globally configures logging.
    '''

    logging.warning(message, *args, **kwargs)


@type_check
def log_error(message: str, *args, **kwargs) -> None:
    '''
    Log the passed error message with the root logger, formatted with the
    passed `%`-style positional and keyword arguments.

    This function expects the :class:`LogConfig` class to have been previously
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
        from betse.util.io.log.logenum import LogLevel
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

        # Generator yielding 2-tuples "(exception, traceback)" for all parent
        # exceptions of this exception *AND* this exception (in that order),
        # where "exception" is each exception and "traceback" is the traceback
        # stored for each exception.
        #
        # Sadly, this list is only gettable via the private
        # traceback._iter_chain() function in older versions of Python. Since
        # this function is unavailable in newer versions of Python, a
        # BETSE-specific compatibility function is called instead.
        exc_parents_generator = stderrs._iter_chain(exception, exc_traceback)

        # Tuple of 2-tuples "(exception, traceback)" in the reverse order
        # yielded by this generator, preserving readability by ensuring that
        # this exception is logged first, the parent exception of this exception
        # (if any) is logged second, and so forth.
        exc_parents = tuple(reversed(tuple(exc_parents_generator)))

        # 0-based index of the last exception in this list.
        exc_parent_last_index = len(exc_parents) - 1

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

        # For each parent exception and that exception's traceback...
        for exc_parent_index, (exc_parent, exc_parent_traceback) in (
            enumerate(exc_parents)):
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

            # If the exception type prefixing the last line of this message is
            # itself prefixed by the expected and hence ignorable
            # fully-qualified name of the subpackage defining BETSE exceptions,
            # truncate this prefix for brevity.
            #
            # Note that the format_exception_only() function guarantees the
            # last line of this message to *ALWAYS* be "the message indicating
            # which exception occurred."
            if exc_message_lines[-1].startswith('betse.exceptions.'):
                exc_message_lines[-1] = exc_message_lines[-1][
                    len('betse.exceptions.'):]

            # Last line of this message. By design, the format_exception_only()
            # function guarantees this line to *ALWAYS* be "the message
            # indicating which exception occurred."
            exc_message_line = exc_message_lines[-1]

            # Append this message to the traceback buffer *BEFORE* appending a
            # truncation of this message to the message buffer.
            exc_full_buffer.write(strs.join(exc_message_lines))
            #print('exception string: '+ exc_message_lines[-1])

            # Split the last line of this message into a non-human-readable
            # exception class and ideally human-readable exception message. If
            # this exception is not "None" *AND* is convertable without raising
            # exceptions into a string, both format_exception_only() and
            # _format_final_exc_line() guarantee this line to be formatted as:
            #     "${exc_class}: ${exc_message}"
            assert types.is_sequence_nonstr_nonempty(exc_message_lines), (
                types.assert_not_sequence_nonstr_nonempty(
                    exc_message_lines, 'Exception message lines'))
            exc_message_match_groups = regexes.get_match_groups_numbered(
                exc_message_line, r'^({})(?:\s*|:\s+(.+))$'.format(
                    identifiers.IDENTIFIER_QUALIFIED_REGEX))

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
                    regexes.remove_substrs(
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
            # this exception. Replace this by a human-readable message.
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

            # If this exception is *NOT* the last, append an explanatory header.
            if exc_parent_index != exc_parent_last_index:
                exc_full_buffer.write(
                    '\n'
                    'The above exception wrapped '
                    'the following originating exception:'
                    '\n\n'
                )

        # Append a random error haiku to the traceback buffer... *BECAUSE*!
        exc_full_buffer.write('\n{}'.format(stderrs.get_haiku_random()))

        # String contents of the traceback buffer.
        exc_full = exc_full_buffer.getvalue()

        # Singleton logging configuration for the current Python process.
        log_config = logconfig.get()

        # If the end user requested that nothing be logged to disk, respect this
        # request by logging tracebacks to the error level and hence stderr.
        # (Avoid printing the synopsis already embedded in these tracebacks.)
        if log_config.file_level >= LogLevel.NONE:
            log_error(exc_full)
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
                # Append a reference to the file being logged to to the
                # synopsis buffer.
                exc_iota_buffer.write(
                    '\n\nFor details, see "{}".'.format(log_config.filename))

                # String contents of the synopsis buffer.
                exc_iota = exc_iota_buffer.getvalue()

                # Print this synopsis.
                stderrs.output(exc_iota)

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
                log_debug(exc_full)
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
