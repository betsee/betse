#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level shell facilities.
'''

# ....................{ IMPORTS                            }....................
import os

# ....................{ GETTERS                            }....................
def get_environment() -> dict:
    '''
    Dictionary mapping the name of each **environment variable** (i.e.,
    external variable imported from the current shell environment for the
    active Python interpreter) to the string value of this variable.

    For safety, this dictionary is guaranteed to be a deep copy of the current
    environment and hence may be safely:

    * Modified in this process _without_ modifying the actual environment.
    * Isolate this environment to subprocesses, preventing concurrent changes in
      the environment of this process from affecting these subprocesses.
    '''

    return os.environ.copy()
