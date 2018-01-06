#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
BETSE-specific test exception hierarchy.
'''

# ....................{ IMPORTS                            }....................
from abc import ABCMeta

# ....................{ EXCEPTIONS                         }....................
class BetseTestException(Exception, metaclass=ABCMeta):
    '''
    Abstract base class of all BETSE-specific test exceptions.
    '''
    pass


class BetseTestFixtureException(BetseTestException):
    '''
    Fixture-specific test exception.
    '''
    pass


class BetseTestHookException(BetseTestException):
    '''
    Hook-specific test exception.
    '''
    pass


class BetseTestParamException(BetseTestException):
    '''
    Parameter-specific test exception.
    '''
    pass
