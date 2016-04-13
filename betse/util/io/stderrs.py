#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level standard error facilities.
'''

# ....................{ IMPORTS                            }....................
import random, sys

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
    Get a random haiku to be printed in the event of fatal errors.
    '''
    return '\n'.join(random.choice(HAIKU))

# ....................{ OUTPUTTERS                         }....................
def output(*objects) -> None:
    '''
    Print all passed objects to standard error *without* logging such objects.

    This function is *not* named `print`, as doing so induces spurious errors
    elsewhere.
    '''
    print(*objects, file = sys.stderr)

# --------------------( WASTELANDS                         )--------------------
# from betse.util.type import strs
# from betse.io.file.log import logger
# ....................{ CLASSES                            }....................
# class LoggingPrinter(object):
#     '''
#     Low-level printing and logging of arbitrary objects.
#
#     Attributes
#     ----------
#     '''
#
#     #FUXME: Cache a logging implementation here as a new private field
#     #"_logger".
#
#     def __init__(self):
#         #FUXME: Actually, this should just be a global variable of the new
#         #"betse.io.file.log" module.
#         #FUXME: Right. And since we no longer need this as a field, there's
#         #really no justification for having this as a class. Revert back to
#         #simple functions, please.

# # ....................{ SINGLETONS                         }....................
# printer = LoggingPrinter()
# '''
# Singleton instance of LoggingPrinter(), simplifying printing and logging of
# arbitrary objects in both the CLI and GUI interfaces.
# '''

# Utility functions manipulating **standard error** (i.e., the canonical file
# handle to which human-readable errors are written).
