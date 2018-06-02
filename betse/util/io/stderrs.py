#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level standard error facilities.
'''

#FIXME: Enable the default fault handler as follows at a suitably early point
#in application startup:
#
#    import faulthandler
#    faulthandler.enable()
#
#The above logic prints a detailed traceback on all segmentation faults to
#standard error, as sanity requires. Note that, while the faulthandler.enable()
#function called above technically supports redirection to an open file handle,
#the handle is required to remain open for the lifetime of the fault handler
#(i.e., foreover). This isn't really feasible, in general, so let's just run
#with the default logic for safety.
#
#Ideally, this submodule would define a new enable_fault_handler() function
#trivially calling faulthandler.enable(); the "betse.ignition" submodule should
#then call that function at a suitably early point -- possibly even as the
#first call, given the horror of segmentation faults that lack tracebacks.
#
#Lastly, note that there are *NO* performance penalties of enabling this
#handler and that effectively all Python applications should do so. For further
#details, see this authoritative StackOverflow answer by the author of the
#"faulthandler" module:
#    https://stackoverflow.com/a/29246977/2809027

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
