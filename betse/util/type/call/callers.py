#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
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
from betse.exceptions import BetseCallableException
from betse.util.type.types import type_check, CallableTypes

# ....................{ GETTERS                            }....................
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

# ....................{ GETTERS ~ caller : basename        }....................
@type_check
def get_caller_basename(call_stack_index: int = 2) -> str:
    '''
    **Basename** (i.e., unqualified name *not* preceded by the ``.``-delimited
    name of the parent module or class) of the callable with the passed index on
    the call stack if any *or* raise an exception otherwise.

    Parameters
    ----------
    call_stack_index : optional[int]
        0-based index of the call stack frame to be inspected. Equivalently, the
        1-based number of leading (most recent) stack frames to be ignored
        including the call to this function. Defaults to an index inspecting the
        **caller's caller** (i.e., the function or method calling the function
        or method calling this function).

    Returns
    ----------
    str
        Basename of this callable, equivalent to the callable after ignoring the
        passed number of leading stack frames on the call stack. As example, if
        passed 2, this is the basename of the second callable on the call stack
        corresponding to that of the caller's caller.

    Raises
    ----------
    BetseCallableException
        If this index exceeds the length of the call stack.
    '''

    # Call stack, including the call to this function.
    call_stack = inspect.stack()

    # Caller stack frame and frame metadata, nullified to avoid exceptions on
    # deleting these locals in the "finally" block below.
    caller_frame = None
    caller_frame_metadata = None

    # Attempt to...
    try:
        # Last valid index into the call stack.
        call_stack_index_max = len(call_stack) - 1

        # If no such callable exists, raise an exception.
        if call_stack_index > call_stack_index_max:
            raise BetseCallableException(
                'Call stack frame {} not found '
                '(i.e., not in range [0, {}]).'.format(
                    call_stack_index, call_stack_index_max))

        # Caller stack frame metadata as an instance of the 5-tuple
        # "(frame, filename, lineno, function, code_context, index)".
        caller_frame_metadata = call_stack[call_stack_index]

        # Caller stack frame as a "FrameType" instance with these fields:
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

        # Return the basename of the callable associated with this frame's code
        # object. Lo, and the "inspect" module's API doth mightily sucketh!
        return caller_frame.f_code.co_name
        # return caller_frame.f_func.__qualname__
    # For safety, explicitly release *ALL* call stack frames obtained above.
    # Failing to do so invites memory leaks due to circular references. See:
    #     https://docs.python.org/3/library/inspect.html#the-interpreter-stack
    finally:
        del call_stack, caller_frame_metadata, caller_frame


@type_check
def get_caller_basename_first_matching(predicate: CallableTypes) -> str:
    '''
    Basename of the first callable whose basename matches the passed predicate
    in the call stack if any *or* raise an exception otherwise.

    Parameters
    ----------
    predicate : CallableTypes
        Callable iteratively passed the basename of each callable on the call
        stack (starting at the call to this function and iterating up the call
        stack), returning ``True`` only if that basename matches this predicate.

    Returns
    ----------
    str
        Basename of the first callable matching this predicate.

    Raises
    ----------
    BetseCallableException
        If no callable on the call stack matches this predicate.

    See Also
    ----------
    :func:`get_caller_basename`
        Further details.
    '''

    # Call stack, including the call to this function.
    call_stack = inspect.stack()

    # Caller stack frame and frame metadata, nullified to avoid exceptions on
    # deleting these locals in the "finally" block below.
    caller_frame = None
    caller_frame_metadata = None

    # Attempt to...
    try:
        # For each stack frame on the call stack as an instance of the 5-tuple
        # "(frame, filename, lineno, function, code_context, index)"...
        for caller_frame_metadata in call_stack:
            # Caller stack frame as a "FrameType" instance with these fields:
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

            # Caller stack frame code object if any or None otherwise.
            caller_frame_code = getattr(caller_frame, 'f_code', None)

            # If this frame is associated with a code object...
            if caller_frame_code is not None:
                # Basename of the callable associated with this frame's code
                # object if any or None otherwise.
                caller_basename = getattr(caller_frame_code, 'co_name', None)

                # If this basename both exists and matches this predicate,
                # return this basename.
                if caller_basename is not None and predicate(caller_basename):
                    return caller_basename
        # Else, no callable on the call stack matches this predicate. In this
        # case, raise an exception.
        else:
            raise BetseCallableException(
                'Caller basename matching predicate {!r} not found.'.format(
                    predicate))
    # For safety, explicitly release *ALL* call stack frames obtained above.
    # Failing to do so invites memory leaks due to circular references. See:
    #     https://docs.python.org/3/library/inspect.html#the-interpreter-stack
    finally:
        del call_stack, caller_frame_metadata, caller_frame

# ....................{ GETTERS ~ caller : module name     }....................
@type_check
def get_caller_module_name(call_stack_index: int = 2) -> str:
    '''
    Fully-qualified name of the parent module of the callable with the passed
    index on the call stack if any *or* raise an exception otherwise.

    Parameters
    ----------
    call_stack_index : optional[int]
        0-based index of the call stack frame to be inspected. Equivalently, the
        1-based number of leading (most recent) stack frames to be ignored
        including the call to this function. Defaults to an index inspecting the
        **caller's caller** (i.e., the function or method calling the function
        or method calling this function).

    Returns
    ----------
    str
        Fully-qualified name of the parent module of this callable, equivalent
        to the callable after ignoring the passed number of leading stack frames
        on the call stack. As example, if passed 2, this is the name of the
        parent module of the second callable on the call stack corresponding to
        that of the caller's caller.

    Raises
    ----------
    BetseCallableException
        If this index exceeds the length of the call stack.

    See Also
    ----------
    https://www.calazan.com/how-to-retrieve-the-name-of-the-calling-module-in-python
        Article strongly inspiring this implementation.
    '''

    # Call stack, including the call to this function.
    call_stack = inspect.stack()

    # Caller stack frame metadata, nullified to avoid exceptions on deleting
    # these locals in the "finally" block below.
    caller_frame_metadata = None

    # Attempt to...
    try:
        # Last valid index into the call stack.
        call_stack_index_max = len(call_stack) - 1

        # If no such callable exists, raise an exception.
        if call_stack_index > call_stack_index_max:
            raise BetseCallableException(
                'Call stack frame {} not found '
                '(i.e., not in range [0, {}]).'.format(
                    call_stack_index, call_stack_index_max))

        # Caller stack frame metadata as an instance of the 5-tuple
        # "(frame, filename, lineno, function, code_context, index)".
        caller_frame_metadata = call_stack[call_stack_index]

        # Return the fully-qualified name of the parent module of the callable
        # associated with this frame.
        return inspect.getmodulename(caller_frame_metadata[1])
    # For safety, explicitly release *ALL* call stack frames obtained above.
    # Failing to do so invites memory leaks due to circular references. See:
    #     https://docs.python.org/3/library/inspect.html#the-interpreter-stack
    finally:
        del call_stack, caller_frame_metadata
