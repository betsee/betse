#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level standard error facilities.
'''

# ....................{ IMPORTS                           }....................
import sys, traceback

# ....................{ OUTPUTTERS                        }....................
def output(*objects) -> None:
    '''
    Print all passed objects to stderr *without* logging these objects.

    This function is intentionally *not* named :func:`print`. Doing so
    introduces subtle issues elsewhere.
    '''

    print(*objects, file=sys.stderr)


def output_exception(heading: str = None) -> None:
    '''
    Print the currently caught exception to stderr *without* logging this
    exception optionally preceded by the passed human-readable heading if any.

    Parameters
    ----------
    heading : optional[str]
        Optional human-readable heading to be printed before this exception if
        any *or* ``None`` if no heading is to be printed.
    '''

    #FIXME: Assert that an exception has actually been raised here.

    # If a heanding is passed, print this heading to stderr.
    if heading is not None:
        output(heading)

    # Print this exception to stderr.
    traceback.print_exc(file=sys.stderr)


def output_traceback() -> None:
    '''
    Print the current call stack to stderr *without* logging this call stack.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.call import callers

    # Print this call stack, excluding the calls to both this and the
    # callers.get_traceback() functions.
    output(callers.get_traceback(-2))
