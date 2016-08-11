#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Unit tests for BETSE's scripting API.
'''

#FIXME: Provide additional tests actually exercising scripting-specific
#functionality exposed by this API (e.g., the betse.script.api.seed() function).
#Doing so requires the customary creation of a test simulation configuration.
#Unfortunately, the fixtures for creating such configurations currently reside
#under the functional test-specific "betse_test.func.fixtures" subpackage rather
#than the more general-purpose "betse_test.fixtures" subpackage. It appears that
#such fixtures (e.g., "betse_sim_config") must now be moved from the former into
#the latter. In theory, this shouldn't be terribly arduous.

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
