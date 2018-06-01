#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
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

# ....................{ IMPORTS                            }....................
import random, sys, traceback

# ....................{ CONSTANTS                          }....................
#FIXME: Shift into a new "betse.util.io.haiku" submodule.
HAIKU = [
    (
        'Chaos reigns within.',
        'Reflect, repent, and reboot.',
        'Order shall return.',
    ), (
        'ABORTED effort:',
        'Close all that you have.',
        'You ask way too much.',
    ), (
        'First snow, then silence.',
        'This thousand dollar screen dies',
        'so beautifully.',
    ), (
        'A crash reduces',
        'your expensive computer',
        'to a simple stone.',
    ), (
        'Error messages',
        'cannot completely convey.',
        'We now know shared loss.',
    ), (
        'The code was willing.',
        'It considered your request',
        'but the chips were weak.',
    ), (
        'There is a chasm',
        'of carbon and silicon',
        "this software can't bridge.",
    ), (
        'To have no errors',
        'would be life without meaning.',
        'No struggle, no joy.',
    ), (
        'many fingers clicking',
        'screens are full of letters',
        'what is their meaning?',
    ), (
        'Water spills downwards.',
        'Electric stream cascades sparks.',
        'Data flow ceases.',
    ), (
        'Technical support',
        'would be a flowing source of',
        'sweet commiseration.',
    ), (
        'Fatal exception.',
        'Code has looped upon itself',
        'like the coiled serpent.',
    ), (
        'An old CPU.',
        'A new program is loaded.',
        'The sound of crashing.',
    ),
]
'''
List of haikus to be printed in the event of fatal errors.

This is serious business, folks.

See Also
----------
http://baetzler.de/humor/haiku_error.var
https://www.gnu.org/fun/jokes/error-haiku.html
    To quote: "IMAGINE IF INSTEAD OF CRYPTIC TEXT STRINGS, YOUR COMPUTER
    PRODUCED ERROR MESSAGES IN HAIKU..." We need no longer imagine.
'''

# ....................{ GETTERS                            }....................
#FIXME: Shift into a new "betse.util.io.haiku" submodule.
def get_haiku_random() -> str:
    '''
    Random instructive haiku to be printed in the event of fatal errors.
    '''

    return '\n'.join(random.choice(HAIKU))


# def get_exception_traceback

# ....................{ OUTPUTTERS                         }....................
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
