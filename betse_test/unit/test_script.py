#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Unit tests for BETSE's scripting API.
'''

#FIXME: Provide additional tests actually exercising scripting-specific
#functionality exposed by this API (e.g., the betse.script.api.seed() function).
#Doing so requires the customary creation of a test simulation configuration,
#available via the general-purpose "betse_sim_conf" fixture.

# ....................{ IMPORTS                            }....................

# ....................{ TESTS                              }....................
def test_script_imports() -> None:
    '''
    Test the importability of BETSE's scripting API.

    This API imports from a medley of other first- and third-party submodules
    and subpackages, which, understandably, occasionally "go AWOL."
    '''

    # Import this API and call all functionality containing these statements.
    from betse import script
