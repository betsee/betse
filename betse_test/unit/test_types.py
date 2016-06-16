#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Unit tests for the `@type_check` decorator, implementing a rudimentary subset of
PEP 484-style type checking based on Python 3.x function annotations.
'''

# ....................{ IMPORTS                            }....................
import pytest

# ....................{ TESTS                              }....................
def test_type_check_noop() -> None:
    '''
    Test type checking for a function with no function annotations, reducing to
    _no_ type checking.
    '''

    # Import this decorator.
    from betse.util.type.types import type_check

    # Type check a function with no function annotations.
    @type_check
    def khorne(gork, mork):
        return gork + mork

    # Call this function and assert the expected return value.
    assert khorne('WAAAGH!', '!HGAAAW') == 'WAAAGH!!HGAAAW'


# ....................{ TESTS ~ pass                       }....................
def test_type_check_pass_keyword_and_positional() -> None:
    '''
    Test type checking for a function passed both annotated positional and
    keyword parameters.
    '''

    # Import this decorator.
    from betse.util.type.types import type_check

    # Type check a function.
    @type_check
    def slaanesh(daemonette: str, keeper_of_secrets: str) -> str:
        return daemonette + keeper_of_secrets

    # Call this function with both positional and keyword arguments and assert
    # the expected return value.
    assert slaanesh(
        'Seeker of Decadence', keeper_of_secrets="N'Kari") == (
        "Seeker of DecadenceN'Kari")


def test_type_check_pass_keyword_only() -> None:
    '''
    Test type checking for a function passed annotated **keyword-only
    parameters** (i.e., parameters following an `*` or `*args` parameter).
    '''

    # Import this decorator.
    from betse.util.type.types import type_check

    # Type check a function.
    @type_check
    def changer_of_ways(sky_shark: str, *, chaos_spawn: str) -> str:
        return sky_shark + chaos_spawn

    # Call this function with keyword arguments and assert the expected return
    # value.
    assert changer_of_ways(
        'Screamers', chaos_spawn="Mith'an'driarkh") == (
        "ScreamersMith'an'driarkh")

# ....................{ TESTS ~ fail                       }....................
def test_type_check_fail_param() -> None:
    '''
    Test type checking for an annotated function call failing a parameter type
    check.
    '''

    # Import this decorator.
    from betse.util.type.types import type_check

    # Type check a function.
    @type_check
    def eldar(isha: str, asuryan: str) -> str:
        return isha + asuryan

    # Call this function with an invalid type and assert the expected exception.
    with pytest.raises(AssertionError):
        eldar('Mother of the Eldar', 100.100)


def test_type_check_fail_return() -> None:
    '''
    Test type checking for an annotated function call failing a return type
    check.
    '''

    # Import this decorator.
    from betse.util.type.types import type_check

    # Type check a function returning an invalid type.
    @type_check
    def necron(star_god: str, old_one: str) -> str:
        return 60e6

    # Call this function and assert the expected exception.
    with pytest.raises(AssertionError):
        necron("C'tan", 'Elder Thing')

# ....................{ TESTS ~ bad                        }....................
def test_type_check_bad_param() -> None:
    '''
    Test type checking for a function with an unsupported parameter annotation.
    '''

    # Import this decorator.
    from betse.util.type.types import type_check

    # Assert the expected exception from attempting to type check a function
    # with a parameter annotation that is *NOT* a type.
    with pytest.raises(TypeError):
        @type_check
        def nurgle(nurgling: str, great_unclean_one: 'Bringer of Poxes') -> str:
            return nurgling + great_unclean_one


def test_type_check_bad_return() -> None:
    '''
    Test type checking for a function with an unsupported return annotation.
    '''

    # Import this decorator.
    from betse.util.type.types import type_check

    # Assert the expected exception from attempting to type check a function
    # with a return annotation that is *NOT* a type.
    with pytest.raises(TypeError):
        @type_check
        def tzeentch(disc: str, lord_of_change: str) -> 'Player of Games':
            return disc + lord_of_change
