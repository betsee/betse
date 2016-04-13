#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
**Call stack** (i.e., list of objects synopsizing all function and method calls
leading to the current function or method call) facilities.

Caveats
----------
Call stack inspection typically assumes the current Python interpreter to be
**CPython** (i.e., the official Python interpreter) rather than a third-party
alternative (e.g., PyPy, Pyston), implying that functions defined by this module
should typically _not_ be called from production code. Call these functions
_only_ where necessary.
'''

# ....................{ IMPORTS                            }....................
import inspect
from betse.exceptions import BetseExceptionFunction
from betse.util.type import types

# ....................{ GETTERS                            }....................
#FIXME: Contribute solution, once working, back to:
#    https://stackoverflow.com/questions/2654113/python-how-to-get-the-callers-method-name-in-the-called-method
#While numerous answers exist, none leverage the Python 3-specific
#"__qualname__" facility -- dramatically simplifying use. Also, prevailing
#answers fail to safely release stack frames.

def get_caller_name(call_stack_index = 2) -> str:
    '''
    Get the fully-qualified name of the desired function or method on the
    current call stack.

    In this case, "desired" means the function or method with the passed index
    in the call stack -- that is, the function or method obtained _after_
    ignoring the passed number of leading stack frames in the call stack.

    Parameters
    ----------------------------
    call_stack_index : int
        0-based index of the call stack frame to be inspected. Equivalently, the
        1-based number of leading (most recent) stack frames to be ignored
        including the call to this function. By default, this index inspects the
        **caller's caller** (i.e., the function or method calling the function
        or method calling this function).
    '''
    assert types.is_int(call_stack_index),\
        types.assert_not_int(call_stack_index)

    # Call stack, including the call to this function.
    call_stack = inspect.stack()

    # If the desired function or method does not exist, raise an exception.
    if call_stack_index >= len(call_stack):
        raise BetseExceptionFunction(
            'Call stack frame {} not found (i.e., not in the range [1, {}]).'.format(
                call_stack_index, len(call_stack)))

    # Caller stack frame.
    caller_frame = call_stack[call_stack_index][0]

    try:
        # Get this fully-qualified name from this frame's function or method
        # object. (Yea, and the "inspect" module API doth sucketh.)
        return caller_frame[3].__qualname__
    # For safety, explicitly release *ALL* call stack frames obtained above.
    # Failing to do so invites memory leaks due to circular references. See:
    #     https://docs.python.org/3/library/inspect.html#the-interpreter-stack
    finally:
        del caller_frame
