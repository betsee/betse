#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level standard output facilities.
'''

#FIXME: Copy setup.util.get_command_output() here after completion.

# ....................{ IMPORTS                           }....................

# ....................{ OUTPUTTERS                        }....................
def output(*objs) -> None:
    '''
    Print all passed objects as is to standard output *without* logging these
    objects.

    This function is provided only for orthogonality with the function of the
    same name defined by the :mod:`betse.util.io.stderrs` module.

    This function is intentionally *not* named ``print()`` to avoid conflict
    with the builtin function of the same name.
    '''

    print(*objs)


def output_lines(*objs) -> None:
    '''
    Print all passed objects as is to standard output *without* logging these
    objects, delimiting each such object by a newline.
    '''

    print('\n'.join(objs))


def output_traceback() -> None:
    '''
    Print the current call stack to standard output.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.call import callers

    # Print this call stack, excluding the calls to both this and the
    # callers.get_traceback() functions.
    output(callers.get_traceback(-2))
