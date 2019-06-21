#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level **generator** (i.e., objects satisfying the standard generator API)
facilities.
'''

# ....................{ IMPORTS                           }....................
from betse.util.type.types import type_check, GeneratorType, IterableTypes

# ....................{ GENERATORS                        }....................
def empty_generator() -> GeneratorType:
    '''
    Empty generator yielding... absolutely nothing.

    See Also
    ----------
    https://stackoverflow.com/a/36658865/2809027
        StackOverflow answer strongly inspiring this implementation.
    '''

    yield from ()

# ....................{ TESTERS                           }....................
def is_generator(obj: object) -> bool:
    '''
    ``True`` only if the passed object is a generator.
    '''

    # Best things in life are free.
    return isinstance(obj, GeneratorType)

# ....................{ GETTERS                           }....................
@type_check
def get_length(generator: GeneratorType) -> int:
    '''
    Length of the passed **finite generator** (i.e., generator guaranteed to
    yield only a finite number of values).

    Caveats
    ----------
    **This function consumes the passed generator,** as a harmful (albeit
    unavoidable) side effect of calculating the length of this generator.

    **This function fails to halt if the passed generator is infinite,** as a
    harmful (albeit unavoidable) side effect of needing to entirely consume
    this generator in order to calculate its length.

    See Also
    ----------
    https://stackoverflow.com/a/31350424
        StackOverflow answer strongly inspiring this implementation.
    '''

    # Convert this generator into a list and return the length of this list.
    #
    # While clearly non-ideal, this approach is well-known to be the most
    # time-efficient solution. The most space-efficient solution is quite
    # different, but also approximately four times slower. Since time is scarce
    # and space is cheap in the general case, we prefer the sane approach.
    return len(list(generator))

# ....................{ CONVERTERS                        }....................
@type_check
def to_tuple_if_generator(iterable: IterableTypes) -> IterableTypes:
    '''
    Tuple of all items iteratively yielded by the passed iterable if this
    iterable is a generator *or* this iterable as is otherwise (i.e., if this
    iterable is *not* a generator).
    '''

    return tuple(iterable) if is_generator(iterable) else iterable
