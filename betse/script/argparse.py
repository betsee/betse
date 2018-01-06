#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

# ....................{ IMPORTS                            }....................
import argparse
from argparse import Namespace
from betse.exceptions import (
    BetseCLIArgException, BetseCLIArgParserException)
from betse.util.type.types import type_check, SequenceTypes

# ....................{ CLASSES                            }....................
#FIXME: This class' current approach to lifecycle management is fragile and,
#frankly, non-Pythonic. The Pythonic approach is to require that an object be:
#
#* Fully and completely initialized on instantiation (e.g., in the __init__()
#  special method).
#* Fully and completely uninitialized by destruction (e.g., by being set to
#  "None" and hence garbage collected).
#
#The Pythonic approach is the simple approach. It works. Sadly, we go to great
#lengths to circumvent that approach here. This class currently permits objects
#to be instantiated in an uninitialized state (which is bad), used in an
#uninitialized state (such as the __len__() special method, which successfully
#returns 0 rather than failing for uninitialized objects -- which is worse), and
#partially but not completely uninitialized (by calling the uninitialize()
#method, which fails to destroy the `self.argv` attribute). Frankly, this is not
#the way.
#
#Consider refactoring this class as follows:
#
#* Requiring tha argument list currently passed to the set_args() method be
#  passed to the __init__() method instead.
#* Replacing all current calls to the uninitialize() method with nullification
#  of the "argv" singleton object defined at the end of this module. There are
#  a few ways to go about this, but the sanest would probably be to define two
#  new module-level functions resembling:
#
#  def set_args(*args: str) -> None:
#      '''
#      Define the current argument list to be passed to BETSE scripts.
#      '''
#
#      globals argv
#      argv = BetseScriptArguments(*args)
#
#
#  def unset_args() -> None:
#      '''
#      Undefine the currently defined argument list if any.
#      '''
#
#      globals argv
#      argv = None
#
#* Removing the:
#  * "is_initialized" boolean.
#  * set_args() method.
#  * uninitialize() method.
#
#After doing so, we could (...and should!) go even further. This entire class
#could simply subclass the "list" type and be replaced with the empty
#implementation: e.g.,
#
#    class BetseScriptArguments(list):
#        pass
#
#Of course, at that point, we'd might as well just replace this class with the
#"list" type everywhere and be done with it. In other words, this was
#overengineered to the hilt. I respect that, however. An overengineer is an
#eager engineer. Praise be to the Great Code Monkey in the Sky, for He shines
#his lurid LCD light on us all!

class BetseScriptArguments(object):
    '''
    Class sanitizing the passing of argument lists to BETSE scripts within an
    arbitrary caller environment.

    Members
    ----------
    argv : list[str]
        List of all argument strings to be passed to BETSE scripts.
    is_initialized : bool
        `True` only if the :meth:`set_args` method has been called, in which
        case the :attr:`argv` list is guaranteed to be non-`None`.
    '''


    def __init__(self) -> None:
        super().__init__()

        self.is_initialized = False
        self.argv = []


    def __len__(self) -> int:
        '''
        Get the size of the previously set argument list if any _or_ 0
        otherwise.

        Returns
        ----------
        int
            This argument list's size.
        '''

        return len(self.argv) if self.is_initialized else 0


    def __getindex__(self, index: int) -> str:
        '''
        Get the previously set argument string with the passed 0-based index
        _or_ raise an exception if this argument list has yet to be set.

        Parameters
        ----------
        index : int
            0-based index of this argument string.

        Returns
        ----------
        str
            This argument string.

        Raises
        ----------
        BetseCLIArgException
            If this argument list has yet to be set.
        '''

        # If an argument list has yet to be set, raise an exception.
        if not self.is_initialized:
            raise BetseCLIArgException('Argument list undefined.')

        # Else, return the desired argument.
        return self.argv[index]


    @type_check
    def set_args(self, *args: str) -> None:
        '''
        Set the argument list to be passed to subsequent script calls.

        Parameters
        ----------
        args : tuple[str]
            Tuple of all such argument strings, internally converted into a
            proper list by this method.
        '''

        self.is_initialized = True
        self.argv = list(args)


    #FIXME: For safety, shouldn't this set "self.argv = None" as well?
    def uninitialize(self) -> None:
        '''
        Uninitialize this object, preventing scripting logic elsewhere from
        leveraging the argument list previously passed to the :meth:`set_args`
        method if any.
        '''

        self.is_initialized = False


class ArgumentParser(argparse.ArgumentParser):
    '''
    Thin subclass providing argument parsing for BETSE scripts.

    This subclass, together with the `argv` global variable, permits BETSE
    scripts to accept arguments both within the BETSE runtime and as standalone
    scripts using exactly the same call.
    '''


    def __init__(self, prog: str = None, **kwargs):
        '''
        Initialize this argument parser.

        Parameters
        ----------
        prog : optional[str]
            Basename of the current executable. Defaults to `None`. If `None`
            _and_ an argument list was passed (via the `argv` global), this
            basename defaults to the first element of that list.

        See the superclass docstring for all remaining keyword arguments.
        '''

        if prog is None and argv.is_initialized:
            super().__init__(prog=argv[0], **kwargs)
        else:
            super().__init__(prog=prog, **kwargs)


    def parse_args(
        self,
        args: SequenceTypes = None,
        namespace: Namespace = None
    ) -> Namespace:
        '''
        Parse all low-level argument strings in the passed argument list into
        their corresponding high-level objects within the passed namespace and,
        for convenience, return the same namespace.

        Parameters
        ----------
        args : optional[SequenceTypes]
            List of all low-level argument strings to be parsed. Defaults to
            `None`. If `None` _and_ an argument list was passed (via the `argv`
            global), this list defaults to all elements of the current `argv`
            global excluding the first element (i.e., the basename of the
            current executable).
        namespace : optional[Namespace]
            Namespace to add attributes for all high-level objects parsed from
            these low-level strings. Defaults to `None`, in which case a new
            `Namespace` instance is instantiated and returned.

        Returns
        ----------
        Namespace
            Passed `namespace` parameter if non-`None` _or_ a new `Namespace`
            instance otherwise.
        '''

        try:
            if args is None and argv.is_initialized:
                return super().parse_args(args=argv[1:], namespace=namespace)
            return super().parse_args(args=args, namespace=namespace)
        except SystemExit as e:
            raise BetseCLIArgParserException(e.code)

# ....................{ SINGLETONS                         }....................
argv = BetseScriptArguments()
'''
Safe `sys.argv` substitute, permitting scripts run within the BETSE runtime to
safely accept a passed argument list _without_ unsafely overwriting the global
`sys.argv` list.
'''
