#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level string facilities.
'''

# ....................{ IMPORTS                            }....................
from betse.util.type import containers

# ....................{ TESTERS                            }....................
def is_str(obj) -> bool:
    '''
    True if the passed object is a *string* (i.e., is an instance of the `str`
    class or a subclass thereof).
    '''
    return isinstance(obj, str)

# ....................{ JOINERS                            }....................
def join(*strings) -> str:
    '''
    Concatenate the passed strings with no separating delimiter.

    This is a convenience function wrapping the standard `"".join((...))`
    method, whose syntax is arguably overly verbose.
    '''
    return join_on(strings, delimiter = '')

def join_on(*strings, delimiter) -> str:
    '''
    Join the passed strings with the passed separating delimiter.

    This is a convenience function wrapping the standard
    `"...".join((...))` method, whose syntax is arguably overly verbose.
    '''
    assert isinstance(delimiter, str), '"{}" not a string.'.format(delimiter)

    # If only one object was passed *AND* such object is a non-string iterable,
    # set the list of passed strings to such object.
    if len(strings) == 1 and containers.is_iterable_nonstring(strings[0]):
        strings = strings[0]

    # Join such strings.
    return delimiter.join(strings)

def join_on_newline(*strings) -> str:
    '''
    Join the passed strings with newline as the separating delimiter.

    This is a convnience function wrapping the standard
    `"\n".join((...))` method, whose syntax is arguably overly verbose.
    '''
    return join_on(*strings, delimiter = '\n')

# --------------------( WASTELANDS                         )--------------------
