#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level exception handling facilities.
'''

# ....................{ IMPORTS                            }....................
import traceback
from betse.util.type import types
from betse.util.type.types import type_check
from io import StringIO

# ....................{ GETTERS                            }....................
@type_check
def get_traceback(exception: Exception) -> str:
    '''
    Non-human-readable traceback associated with the passed exception, typically
    intended to be logged and/or displayed to end users for debugging purposes.

    Parameters
    ----------
    exception : Exception
        Exception to return this traceback for.

    Returns
    ----------
    str
        Non-human-readable traceback associated with this exception.
    '''

    # Terse synopsis and verbose traceback for this exception.
    _, exc_traceback = get_metadata(exception)

    # Return this traceback.
    return exc_traceback


@type_check
def get_metadata(exception: Exception) -> tuple:
    '''
    Tuple of various metadata specific to the passed exception, typically
    intended to be logged and/or displayed to end users.

    Parameters
    ----------
    exception : Exception
        Exception to return metadata for.

    Returns
    ----------
    (str, str)
        2-tuple ``(synopsis, traceback)``, where:
        * ``synopsis`` is this exception's human-readable synopsis.
        * ``traceback`` is this exception's non-human-readable traceback.
    '''

    # Avoid circular import dependencies.
    from betse.util.io import stderrs
    from betse.util.py import pyident
    from betse.util.type.text import regexes, strs

    # Generator yielding 2-tuples "(exception, traceback)" for all parent
    # exceptions of this exception *AND* this exception (in that order), where
    # "exception" is each exception and "traceback" is the traceback stored for
    # each exception.
    #
    # Sadly, this list is only gettable via the private traceback._iter_chain()
    # function in older versions of Python. Since this function is unavailable
    # in newer versions of Python, an application-specific compatibility
    # function is called instead.
    exc_parents_generator = _iter_chain(exception, exception.__traceback__)

    # Tuple of 2-tuples "(exception, traceback)" in the reverse order yielded by
    # this generator, preserving readability by ensuring that this exception is
    # logged first, the parent exception of this exception (if any) is logged
    # second, and so forth.
    exc_parents = tuple(reversed(tuple(exc_parents_generator)))

    # 0-based index of the last exception in this list.
    exc_parent_last_index = len(exc_parents) - 1

    # String buffer containing a human-readable synopsis of each exception in
    # this chain, unconditionally output to stderr.
    exc_iota_buffer = StringIO()

    # String buffer containing a non-human-readable traceback of each exception
    # in this chain, conditionally logged to the logfile.
    exc_full_buffer = StringIO()

    # Human-readable header prefixing each such buffer.
    buffer_header = 'Exiting prematurely due to fatal error:\n\n'

    # Initialize these buffers to this header.
    exc_iota_buffer.write(buffer_header)
    exc_full_buffer.write(buffer_header)

    # For each parent exception and that exception's traceback...
    for exc_parent_index, (exc_parent, exc_parent_traceback) in (
        enumerate(exc_parents)):
        # If this exception is a string, append this string to the synopsis
        # buffer as is and continue to the next parent. This is an edge case
        # that should *NEVER* happen... but could.
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
        assert types.is_sequence_nonstr_nonempty(exc_message_lines), (
            types.assert_not_sequence_nonstr_nonempty(
                exc_message_lines, 'Exception message lines'))

        # If the exception type prefixing the last line of this message is
        # itself prefixed by the expected and hence ignorable fully-qualified
        # name of the subpackage defining BETSE exceptions, truncate this prefix
        # for brevity.
        #
        # Note that the format_exception_only() function guarantees the last
        # line of this message to *ALWAYS* be "the message indicating which
        # exception occurred."
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
        # exception class and ideally human-readable exception message. If this
        # exception is not "None" *AND* is convertable without raising
        # exceptions into a string, both format_exception_only() and
        # _format_final_exc_line() guarantee this line to be formatted as:
        #     "${exc_class}: ${exc_message}"
        exc_message_match_groups = regexes.get_match_groups_numbered(
            exc_message_line, r'^({})(?:\s*|:\s+(.+))$'.format(
                pyident.IDENTIFIER_QUALIFIED_REGEX))

        # This message is guaranteed to be prefixed by a class name.
        exc_class_name = exc_message_match_groups[0]

        # This message is *NOT* guaranteed to be prefixed by a non-empty message
        # (e.g., assert statements passed no message).
        #
        # If a non-empty message matched, use that.
        exc_message = None
        if exc_message_match_groups[1] is not None:
            exc_message = exc_message_match_groups[1]
        # Else if a debug assertion failed with no explicit message, use the
        # exception context directly detailing this assertion.
        elif exc_class_name == 'AssertionError':
            exc_message = 'Debug assertion failed: {}'.format(
                # A traceback line typically contains an internal newline. The
                # substring preceding this newline details the file and function
                # containing the corresponding call; the substring following
                # this newline is this call. Hence, ignore the former.
                regexes.remove_substrs(
                    exc_traceback_lines[-1], r'^.+\n\s*'))
        # Else, convert this exception's class name into a human-readable
        # message (e.g., from "FileNotFoundError" to "File not found error.").
        # Well, try... at least!
        else:
            exc_message = strs.uppercase_char_first(
                pyident.convert_camelcase_to_whitespaced_lowercase(
                    exc_class_name))
        assert types.is_str_nonempty(exc_message), (
            types.assert_not_str_nonempty(
                exc_message, 'Exception message'))

        # If this class is "KeyError", this message is the single-quoted name of
        # a non-existent key in a dictionary whose access raised this exception.
        # Replace this by a human-readable message.
        if exc_class_name == 'KeyError':
            exc_message = 'Dictionary key {} not found.'.format(
                exc_message)

        # Append this message to the synopsis buffer. For readability, this
        # message is wrapped to the default terminal width and each wrapped line
        # prefixed by indentation.
        exc_iota_buffer.write(strs.wrap(
            text=exc_message, line_prefix='    '))

        # If this exception has a traceback, append this traceback to the
        # traceback but *NOT* synopsis buffer.
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

    # Return the string contents of these buffers in the expected order.
    return (exc_iota_buffer.getvalue(), exc_full_buffer.getvalue())

# ....................{ PRIVATE ~ iterators                }....................
# If the active Python interpreter is 3.4, import the private _iter_chain()
# method from the standard "traceback" module.
try:
    from traceback import _iter_chain
    if False: _iter_chain  # squelch IDE warnings
# Else, the active Python interpreter is >= 3.5, which replaced this method with
# a new public class hierarchy (e.g., "TracebackException"). For portability,
# forward port the traceback._iter_chain() method from the most recent stable
# release of Python 3.4.
#
# Alternately, this class hierarchy *COULD* be backported from the most recent
# stable release of Python 3.6. Doing so, however, would be considerably more
# difficult than simply defining a single function. Thus the current approach.
except ImportError:
    # Private string constants required by _iter_chain(). Fortuitously, these
    # constants remain unchanged in Python >= 3.5.
    from traceback import _cause_message, _context_message

    def _iter_chain(exc, custom_tb=None, seen=None):
        '''
        Private :func:`traceback._iter_chain` method forward-ported without
        modification from the most recent stable release of Python 3.4 as of
        this writing (i.e., Python 3.4.3).
        '''

        if seen is None:
            seen = set()
        seen.add(exc)
        its = []
        context = exc.__context__
        cause = exc.__cause__
        if cause is not None and cause not in seen:
            its.append(_iter_chain(cause, False, seen))
            its.append([(_cause_message, None)])
        elif (context is not None and
            not exc.__suppress_context__ and
            context not in seen):
            its.append(_iter_chain(context, None, seen))
            its.append([(_context_message, None)])
        its.append([(exc, custom_tb or exc.__traceback__)])
        # itertools.chain is in an extension module and may be unavailable
        for it in its:
            yield from it
