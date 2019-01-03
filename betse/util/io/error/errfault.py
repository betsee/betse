#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
**Segmentation fault** (i.e., memory access violation by the active Python
process, caused by attempting to access an inaccessible, illegal, or otherwise
invalid memory location) facilities.
'''

# ....................{ IMPORTS                           }....................
import faulthandler

# ....................{ HANDLERS                          }....................
def handle_faults() -> None:
    '''
    Enable Python's standard handler for **segmentation faults** (i.e., memory
    access violations by the active Python process, caused by attempting to
    access inaccessible, illegal, or otherwise invalid memory locations).

    On the first segmentation fault, this handler prints a detailed traceback
    of the current call stack for the current thread to standard error before
    abruptly terminating this process.

    Caveats
    ----------
    Technically, the :func:`faulthandler.enable` function internally called by
    this function supports redirection to an open file handle (e.g., referring
    to a user-specific log file) rather than to standard error. However, that
    handle is required to remain open for the lifetime of the fault handler and
    hence the active Python process. Since guaranteeing this is infeasible in
    practice, this function defers to the default fault handler for safety.

    Note that there are *no* performance penalties associated with enabling
    this handler. Ergo, *all* Python applications are advised to do so. For
    further details, see `this authoritative StackOverflow answer`_ by the
    author of the standard :mod:`faulthandler` module.

    .. _authoritative StackOverflow answer:
       https://stackoverflow.com/a/29246977/2809027

    See Also
    ----------
    :func:`ignore_faults`
        Function disabling this handler, thus reverting to Python's default
        behavior with respect to segmentation faults.
    '''

    faulthandler.enable()


def ignore_faults() -> None:
    '''
    Disable Python's standard handler for **segmentation faults** (i.e., memory
    access violations by the active Python process, caused by attempting to
    access inaccessible, illegal, or otherwise invalid memory locations).

    This function reverts to Python's default behavior with respect to
    segmentation faults. On the first segmentation fault, Python will abruptly
    terminate this process with *no* traceback or advance notice.

    See Also
    ----------
    :func:`handle_faults`
        Function re-enabling this handler.
    '''

    faulthandler.disable()
