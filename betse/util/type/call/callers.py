#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level **call stack** (i.e., list of objects synopsizing all function and
method calls leading to the current function or method call) functionality.

Caveats
----------
Call stack inspection typically assumes the current Python interpreter to be
**CPython** (i.e., the official Python interpreter) rather than a third-party
alternative (e.g., PyPy, Pyston), implying that functions defined by this module
should typically *never* be called by production code. Consider calling these
functions only where needed.

See Also
----------
:func:`betse.util.io.stderr.output_traceback`
    Printing the current call stack to standard error.
:func:`betse.util.io.stdout.output_traceback`
    Printing the current call stack to standard output.
'''

# ....................{ IMPORTS                            }....................
import inspect, traceback
from betse.exceptions import BetseFunctionException
from betse.util.type.types import type_check

# ....................{ GETTERS                            }....................
@type_check
def get_caller_basename(call_stack_index: int = 2) -> str:
    '''
    Basename of the callable with the passed index in the call stack.

    This function returns the basename of the callableo obtained *after*
    ignoring the passed number of leading stack frames on the call stack.

    Parameters
    ----------
    call_stack_index : optional[int]
        0-based index of the call stack frame to be inspected. Equivalently, the
        1-based number of leading (most recent) stack frames to be ignored
        including the call to this function. Defaults to an index inspecting the
        **caller's caller** (i.e., the function or method calling the function
        or method calling this function).
    '''

    # Call stack, including the call to this function.
    call_stack = inspect.stack()

    # If the desired function or method does not exist, raise an exception.
    if call_stack_index >= len(call_stack):
        raise BetseFunctionException(
            'Call stack frame {} not found '
            '(i.e., not in the range [1, {}]).'.format(
                call_stack_index, len(call_stack)))

    # Caller stack frame metadata as an instance of the 5-tuple
    # "(frame, filename, lineno, function, code_context, index)".
    caller_frame_metadata = call_stack[call_stack_index]

    # Caller stack frame as a "FrameType" instance providing these fields:
    #
    # * "f_back", next outer frame object (this frame's caller).
    # * "f_builtins", built-in namespace seen by this frame.
    # * "f_code", code object being executed in this frame.
    # * "f_globals", global namespace seen by this frame.
    # * "f_lasti", index of last attempted instruction in bytecode.
    # * "f_lineno", current line number in Python source code.
    # * "f_locals", local namespace seen by this frame.
    # * "f_trace", tracing function for this frame or "None".
    caller_frame = caller_frame_metadata[0]

    try:
        # Return the fully-qualified name of this frame's function or method
        # object. Lo, and the "inspect" module's API doth mightily sucketh!
        return caller_frame.f_code.co_name
        # return caller_frame.f_func.__qualname__
    # For safety, explicitly release *ALL* call stack frames obtained above.
    # Failing to do so invites memory leaks due to circular references. See:
    #     https://docs.python.org/3/library/inspect.html#the-interpreter-stack
    finally:
        del call_stack, caller_frame


@type_check
def get_traceback(call_stack_index_last: int = -1) -> str:
    '''
    Human-readable traceback of the current call stack, ignoring exceptions and
    excluding _all_ stack frames inclusively following the stack frame with the
    passed index.

    Parameters
    ----------
    call_stack_index_last : int
        0-based index of the first call stack frame to be ignored. Defaults to
        -1, in which case only the call to this function is ignored.
    '''

    # Call stack excluding the call to the current function.
    call_stack = traceback.extract_stack()[:call_stack_index_last]

    # Human-readable header preceding this traceback.
    traceback_header = 'Traceback (most recent call last):\n'

    # Return this traceback.
    return traceback_header + ''.join(traceback.format_list(call_stack))
