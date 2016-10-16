#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Metadata describing options accepted by BETSE's command line interface (CLI).
'''

# ....................{ IMPORTS                            }....................

# ....................{ OPTIONS                            }....................
#FIXME: Refactor these string globals into a single dictionary mapping from
#option name to help string. The current approach is *MUCH* too heavyweight.

OPTION_VERSION = '''
print program version and exit
'''
'''
Help string template synopsizing the `--version` option.
'''

OPTION_VERBOSE = '''
print low-level debugging messages
'''
'''
Help string template synopsizing the `--verbose` option.
'''

OPTION_LOG_TYPE = '''
type of logging to perform (defaults to "{default}"):
;* "none", logging to stdout and stderr only
;* "file", logging to stdout, stderr, and "--log-file"
'''
'''
Help string template synopsizing the `--log-type` option.
'''

OPTION_LOG_FILE = '''
file to log to if "--log-type" is "file" (defaults to "{default}")
'''
'''
Help string template synopsizing the `--log-file` option.
'''

OPTION_PROFILE_TYPE = '''
type of profiling to perform (defaults to "{default}"):
;* "none", disabling profiling
;* "call", profiling callables (e.g., functions, methods)
;* "line", profiling lines via the third-party "
'''
'''
Help string template synopsizing the `--profile-type` option.
'''

OPTION_PROFILE_FILE = '''
file to profile to if "--profile-type" is not "none" (defaults to "{default}")
'''
'''
Help string template synopsizing the `--profile-file` option.
'''
