#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Unit tests exercising third-party dependencies.
'''

# ....................{ IMPORTS                            }....................
# import pytest

# ....................{ TESTS                              }....................
def test_lib_setuptools_mappings() -> None:
    '''
    Unit test validating each mandatory and optional run- and test-time
    dependency to have a corresponding key of the
    :data:`betse.lib.setuptools.setuptool.SETUPTOOLS_TO_MODULE_NAME` dictionary.
    '''

    # Defer heavyweight imports.
    from betse.metadata import (
        DEPENDENCIES_RUNTIME_MANDATORY,
        DEPENDENCIES_RUNTIME_OPTIONAL,
        DEPENDENCIES_TESTING_MANDATORY,
    )
    from betse.lib.setuptools.setuptool import SETUPTOOLS_TO_MODULE_NAME
    from betse.util.type import mappings

    # Dictionary merged from these dictionaries.
    setuptools_name_to_specs = mappings.merge(
        DEPENDENCIES_RUNTIME_MANDATORY,
        DEPENDENCIES_RUNTIME_OPTIONAL,
        DEPENDENCIES_TESTING_MANDATORY,
    )

    # Assert the setuptools-specific project name of each such dependency to be
    # registered with the "setuptool" submodule.
    for setuptools_name in setuptools_name_to_specs.keys():
        assert setuptools_name in SETUPTOOLS_TO_MODULE_NAME
