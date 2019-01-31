#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level standard error facilities.
'''

# ....................{ IMPORTS                           }....................
import sys, traceback
from betse.util.type.types import type_check, StrOrNoneTypes

# ....................{ OUTPUTTERS                        }....................
@type_check
def output(*texts: str) -> None:
    '''
    Print all passed strings to stderr *without* logging these strings.

    This function is intentionally *not* named :func:`print`, as doing so would
    invite conflicts with the standard :func:`print` function.
    '''

    print(*texts, file=sys.stderr)


@type_check
def output_warning(*warnings: str) -> None:
    '''
    Print all passed strings (prefixed by a suitable human-readable warning
    label) to stderr *without* logging these strings.
    '''

    # ...that was easier than expected.
    output('WARNING: ', *warnings)

# ....................{ OUTPUTTERS ~ exception            }....................
@type_check
def output_exception(heading: StrOrNoneTypes = None) -> None:
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
