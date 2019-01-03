#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level **queue** (i.e., First-In, First-Out (FIFO) data structure, typically
backed by the standard :class:`collections.deque` type) functionality.
'''

# ....................{ IMPORTS                           }....................
import itertools
from betse.util.type.types import type_check, QueueType

# ....................{ POPPERS                           }....................
@type_check
def pop_left(queue: QueueType, count: int) -> None:
    '''
    Pop (i.e., remove) the passed number of items from the head of the passed
    queue in a modestly efficient manner in-place (i.e., modifying the passed
    queue itself).

    Parameters
    ----------
    queue : QueueType
        Queue to pop this number of items from the head of.
    count : int
        Number of items to pop from the head of this queue.

    See Also
    ----------
    https://stackoverflow.com/questions/9507636/how-can-i-pop-lots-of-elements-from-a-deque#comment12041956_9507840
        StackOverflow comment strongly inspiring this implementation.
    '''

    # The following implementation is assumed but has yet to be profiled to be
    # the optimally efficient solution. Alternatives include a map()-based
    # solution, which commonly profiles slower than comparable starmap()-based
    # solutions in the general case:
    #     map(apply, repeat(queue.popleft, count))
    list(itertools.starmap(queue.popleft, itertools.repeat((), count)))
