#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Environmental context for the BETSE REPL.

Each function and variable in this module is loaded into the global
:data:`repl_env` dictionary via a call to :func:`locals`. For safety, this
dictionary is the only attribute that should be imported from this module.
'''

# ....................{ IMPORTS                            }....................
# from betse.util.type.types import type_check

# ....................{ GLOBALS                            }....................
__betse_repl__ = True
'''
`True` only if the current environment is an active REPL.

It is important that scripts be able to tell that they are running in a REPL
as opposed to standalone mode. For example, if the script's author opts to
follow the "if main" convention, then the script will only run in standalone
mode as `__name__ == "betse.repl.environment"` within the REPL. This boolean
simplifies matters by allowing logic resembling:

    if __name__ == '__main__' or __betse_repl__:
        ...

or if behavior should differ between hosted and standalone modes:

    if __name__ == '__main__':
        ...
    elif __betse_repl__:
        ...
'''

# ....................{ GLOBALS ~ env                      }....................
REPL_ENV = locals()
'''
Dictionary mapping from the names to values of all attributes defined by this
submodule, including those defined by the :mod:`betse.script` subpackage.
'''
