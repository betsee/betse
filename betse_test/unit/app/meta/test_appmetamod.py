#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Unit tests exercising functionality defined by the
:mod:`betse.util.app.meta.appmetamod` submodule.
'''

# ....................{ IMPORTS                           }....................
# import pytest
from betse.metadata import VERSION
from betse.util.test.pytest.mark.pytskip import skip_unless_requirement

# ....................{ TESTS                             }....................
# Exercising this function requires at least two modules satisfying the
# "betse.metadata" API. Since the only other module currently doing so is the
# "betsee.guimetadata" submodule, this test requires a recent BETSEE version
# guaranteed to be compatible with the current BETSE version. Since this exact
# merger is the principle use case for this function, exercising this use case
# substantially improves the reliability of BETSE <-> BETSEE integration.
@skip_unless_requirement('BETSEE >= {}.0'.format(VERSION))
def test_merge_module_metadeps() -> None:
    '''
    Unit test the :func:`betse.util.app.meta.appmetamod.merge_module_metadeps`
    function.
    '''

    # Defer heavyweight imports.
    from betse import metadeps as betse_metadeps
    from betse.util.app.meta import appmetamod
    from betse.util.app.meta.appmetamod import MERGE_MODULE_METADEPS_DICTS_NAME
    from betse.util.type.iterable.mapping import maptest
    from betse.util.type.obj import objtest
    from betsee import guimetadeps as betsee_metadeps

    # Application dependency metadata module merging the BETSE and BETSEE
    # application dependency metadata modules.
    module_metadeps_merged = appmetamod.merge_module_metadeps(
        # Fully-qualified name of the module to be created. For simplicity, we
        # create a guaranteeably unique test-specific submodule of the
        # top-level "betse_test" package, which is guaranteed to exist.
        module_name='betse_test.__betse_plus_betsee_metadeps__',
        modules_metadeps=(betse_metadeps, betsee_metadeps,),
    )

    # Test whether this module defines all requisite dictionary globals.
    objtest.die_unless_has_attr(
        module_metadeps_merged, *MERGE_MODULE_METADEPS_DICTS_NAME)

    # Test whether this module defines the following mandatory runtime
    # dependencies:
    #
    # * "Numpy", guaranteed to be required by BETSE.
    # * "PySide2", guaranteed to be required by BETSEE.
    maptest.die_unless_has_keys(
        mapping=module_metadeps_merged.RUNTIME_MANDATORY,
        keys=('Numpy', 'PySide2',))

    # Test whether this module defines the following optional runtime
    # dependencies:
    #
    # * "networkx", guaranteed to be required by BETSE.
    # * "pyside2uic", guaranteed to be required by BETSEE.
    maptest.die_unless_has_keys(
        mapping=module_metadeps_merged.RUNTIME_OPTIONAL,
        keys=('networkx', 'pyside2uic',))

    # Test whether this module defines the following mandatory test-time
    # dependencies:
    #
    # * "pytest", guaranteed to be required by both BETSE and BETSEE.
    maptest.die_unless_has_keys(
        mapping=module_metadeps_merged.TESTING_MANDATORY,
        keys=('pytest',))
