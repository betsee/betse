#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

import argparse
from betse.exceptions import BetseArgumentParserException

class BetseScriptArguments(object):
    '''
    A class to facilitate the passing of arguments to BETSE scripts within a
    calling environment.
    '''
    def __init__(self):
        self.is_initialized = False
        self.argv = []
        super().__init__()

    def __len__(self):
        return len(self.argv) if self.is_initialized else 0

    def __getitem__(self, key):
        return self.argv[key]

    def set_args(self, *args):
        '''
        Set the arguments to be passed to subsequent script calls.
        '''
        self.is_initialized = True
        self.argv = list(args)

    def uninitialize(self):
        '''
        Uninitialize the object so subsequent commands to not attempt to use
        the old arguments.
        '''
        self.is_initialized = False

betse_argv = BetseScriptArguments()
'''
A horrendous global variable to replace sys.argv when calling scripts within
the BETSE runtime.
'''

class ArgumentParser(argparse.ArgumentParser):
    '''
    A thin derived class to provide argument parsing for BETSE scripts.

    This class, together with the `betse_argv` global variable, makes is
    possible for BETSE scripts to accept arguments both within the BETSE
    runtime and as standalone scripts using exactly the same call.
    '''
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def parse_args(self, args=None, namespace=None):
        '''
        Parse the arguments *args* within the *namespace*, returing
        *namespace*.
        '''
        try:
            if args is None and betse_argv.is_initialized:
                return super().parse_args(args=betse_argv, namespace=namespace)
            return super().parse_args(args=args, namespace=namespace)
        except SystemExit as e:
            raise BetseArgumentParserException()
