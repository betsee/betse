#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level standard error facilities.
'''

# ....................{ IMPORTS                            }....................
import random, sys, traceback

# ....................{ CONSTANTS                          }....................
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
def get_haiku_random() -> str:
    '''
    Random haiku to be printed in the event of fatal errors.
    '''

    return '\n'.join(random.choice(HAIKU))

# ....................{ OUTPUTTERS                         }....................
def output(*objects) -> None:
    '''
    Print all passed objects to stderr *without* logging these objects.

    This function is intentionally *not* named `print()`. Doing so introduces
    subtle issues elsewhere.
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
        any _or_ `None` if no heading is to be printed.
    '''

    #FIXME: Assert that an exception has actually been raised here.

    # If a heanding is passed, print this heading to stderr.
    if heading is not None:
        output(heading)

    # Print this exception to stderr.
    traceback.print_exc(file=sys.stderr)


def output_traceback() -> None:
    '''
    Print the current call stack to stderr _without_ logging this call stack.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.call import callers

    # Print this call stack, excluding the calls to both this and the
    # callers.get_traceback() functions.
    output(callers.get_traceback(-2))

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
        Private `traceback._iter_chain()` method forward-ported without
        modification from the most recent stable release of Python 3.4 as of
        this writing: Python 3.4.3.
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
