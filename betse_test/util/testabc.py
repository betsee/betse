#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Abstract base classes (ABCs) globally applicable to all tests.
'''

# ....................{ IMPORTS                            }....................
from abc import ABCMeta  #, abstractmethod

# ....................{ CLASSES                            }....................
class SerialTestABC(metaclass=ABCMeta):
    '''
    Abstract base class running all test methods defined by this concrete
    subclass **serially** (i.e., in the order in which this subclass declares
    these test such that subsequently declared tests depend on the success of
    all previously declared tests).

    On the first failure of such a test, all subsequent such tests will be
    automaticaly marked as xfailing (i.e., failing *without* being run).

    The majority of the black magic required by this class is implemented as
    low-level py.test hooks in the top-level :mod:`betse_func.conftest` plugin.

    FIXME: The following attribute is dynamically set on the py.test object
    encapsulating this class instance rather than this actual class instance.

    Attributes
    ----------
    _first_failure_method_name : str
        Unqualified name of the first failing test method (i.e., method whose
        execution raised an exception) declared by this subclass for the current
        test session if any such method failed _or_ `None` otherwise (i.e., if
        no such methods have yet to fail).
    '''

    # ..................{ INITIALIZERS                       }..................
    @staticmethod
    def is_test_serial(item: 'pytest.main.Item') -> bool:
        '''
        ``True`` only if the passed test callable is **serial** (i.e., a method
        of a subclass of this class).

        Serial methods are intended to be run in test method declaration order,
        such that subsequently declared test methods depend on the success of
        all previously declared test methods.

        Parameters
        ----------
        item : pytest.main.Item
            Metadata encapsulating this test callable (e.g., function, method).

        Returns
        ----------
        bool
            ``True`` only if this test is serial.
        '''

        # Class of this test callable if this callable is a method or "None".
        test_class = item.parent.cls

        # This callable is intended to be run serially only if the test subclass
        # to which this callable is bound subclasses this abstract base class.
        return test_class is not None and issubclass(test_class, SerialTestABC)
