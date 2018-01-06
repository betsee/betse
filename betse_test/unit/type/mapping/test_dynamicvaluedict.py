#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Unit tests exercising the :class:`betse.util.type.mappings.DynamicValueType`
type.
'''

# ....................{ IMPORTS                            }....................
import pytest
import sys
from pytest import fixture

# ....................{ GLOBALS                            }....................
# Global variables whose values are referenced by the dictionary returned by the
# "betse_dynamicvaluedict" fixture.

module = sys.modules[__name__]
'''
Current module object, for use in setattr() calls below.
'''


varda = (
    'O stars that in the Sunless Year',
    'With shining hand by her were sown,',
    'In Windy fields now bright and clear',
    'We see your silver blossom blown!'
)

yavanna = (
    "'All have their worth,' said Yavanna,",
    "'and each contributes to the worth of the others.",
    "But the kelvar can flee or defend themselves,"
    "whereas the olvar that grow cannot.",
    "And among these I hold trees dear.",
    "Long in the growing, swift shall they be in the felling,",
    "and unless they pay toll with fruit upon bough",
    "little mourned in their passing.",
    "So I see in my thought.",
    "Would that the trees might speak on behalf of all things that have roots,",
    "and punish those that wrong them!'",
)

# ....................{ FIXTURES                           }....................
@fixture(scope='session')
def betse_dynamicvaluedict() -> 'DynamicValueDict':
    '''
    Fixture creating and returning a mock :class:`DynamicValueDict` to be
    tested.
    '''

    # Defer heavyweight imports.
    from betse.util.type.mapping.mapcls import DynamicValue, DynamicValueDict

    # Create and return an instance of this type.
    return DynamicValueDict({
        'Varda Elentári': DynamicValue(
            lambda: varda,

            # Sadly, lambda expressions prohibit assignment in Python. *URGH!*
            lambda value: setattr(module, 'varda', value),
            # lambda value: sim.cc_cells.__setindex__(sim.iNa, value)
        ),
        'Yavanna': DynamicValue(
            lambda: yavanna,
            lambda value: setattr(module, 'yavanna', value),
        ),
    })

# ....................{ TESTS                              }....................
def test_dynamicvaluedict_pass(betse_dynamicvaluedict) -> None:
    '''
    Test all aspects of the :class:`DynamicValueDict` type intended to succeed.

    Parameters
    ----------
    betse_dynamicvaluedict : DynamicValueDict
        Instance of this type to be tested.
    '''

    # Test this type's __len__() implementation.
    assert len(betse_dynamicvaluedict) == 2

    # Test this type's __iter__() implementation.
    for dynamicvalue in betse_dynamicvaluedict:
        assert isinstance(dynamicvalue, str)

    # Test this type's __getindex__() implementation.
    assert len(betse_dynamicvaluedict['Varda Elentári']) == 4

    # Test this type's __setindex__() implementation. Specifically, validate
    # that setting a dictionary value sets the value of the global variable
    # underlying this value.
    assert len(yavanna) == 10
    betse_dynamicvaluedict['Yavanna'] = (
        "'Yet it was in the Song,' said Yavanna.",
        "'For while thou wert in the heavens",
        "and with Ulmo built the clouds and poured out the rains,",
        "I lifted up the branches of great trees to receive them,",
        "and some sang to Ilúvatar amid the wind and the rain.'",
    )
    assert len(yavanna) == 5


def test_dynamicvaluedict_fail(betse_dynamicvaluedict) -> None:
    '''
    Test all aspects of the :class:`DynamicValueDict` type intended to fail.

    Parameters
    ----------
    betse_dynamicvaluedict : DynamicValueDict
        Instance of this type to be tested.
    '''

    # Imports deferred for safety.
    from betse.exceptions import BetseMethodUnimplementedException

    # Test this type's __setindex__() implementation. Specifically, validate
    # that setting a dictionary value *NOT* predefined at dictionary
    # initialization time fails.
    with pytest.raises(KeyError):
        betse_dynamicvaluedict['Estë'] = (
            'Grey is her raiment,',
            'and rest her gift.',
        )

    # Test this type's __delindex__() implementation.
    with pytest.raises(BetseMethodUnimplementedException):
        del betse_dynamicvaluedict['Varda Elentári']
