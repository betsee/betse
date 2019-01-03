#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
**Error haiku** (i.e., pseudo-random human-readable general-purpose error
message formatted as a haiku) facilities.
'''

# ....................{ IMPORTS                           }....................
import random, sys, traceback

# ....................{ CONSTANTS                         }....................
ERROR_HAIKU = (
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
)
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

# ....................{ GETTERS                           }....................
def get_random() -> str:
    '''
    Pseudo-random human-readable general-purpose error message formatted as a
    haiku, commonly printed and/or logged in the event of fatal errors.
    '''

    return '\n'.join(random.choice(ERROR_HAIKU))
