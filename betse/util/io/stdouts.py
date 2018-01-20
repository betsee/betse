#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level standard output facilities.
'''

#FIXME: Copy setup.util.get_command_output() here after completion.

# ....................{ IMPORTS                            }....................

# ....................{ OUTPUTTERS                         }....................
def output(*objects) -> None:
    '''
    Print all passed objects to stdout _without_ logging these objects.

    This function is provided only for orthogonality with the function of the
    same name defined by the `stderr` module.

    This function is intentionally _not_ named `print()`. Doing so introduces
    subtle issues elsewhere.
    '''

    print(*objects)


def output_lines(*objects) -> None:
    '''
    Print all passed objects to standard output *without* logging such objects,
    delimiting each such object by a newline.
    '''

    print('\n'.join(*objects))


def output_traceback() -> None:
    '''
    Print the current call stack to standard output.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.call import callers

    # Print this call stack, excluding the calls to both this and the
    # callers.get_traceback() functions.
    output(callers.get_traceback(-2))
