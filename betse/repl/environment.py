#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.
'''
This module provides the environmental context for BETSE REPLs.

Each function and variable in this module is loaded into the `repl_env`
dictionary via a call to `locals`. This is the only symbol that should
be imported from this module.
'''
from betse.script import *
from betse.script.runner import run_script

__betse_repl__ = True

def quit():
    '''
    Gracefully exit the REPL, returning control the the caller.
    '''
    raise SystemExit(0)

repl_env = locals()
