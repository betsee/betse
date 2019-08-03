#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
:mod:`betse.science` initialization functionality.

Motivation
----------
This submodule guarantees that, on the first importation of the
:mod:`betse.science` subpackage, both the current application *and* all
mandatory third-party dependencies of this application have been initialized
with sane defaults.

This submodule is effectively syntactic sugar. While technically unnecessary,
this submodule streamlines interactive use (e.g., from web-based Jupyter
notebooks or CLI-based iPython consoles) by implicitly guaranteeing this
application to be fully usable *without* manual intervention by end users.

This submodule silently reduces to a noop when this application has already
been initialized, as is the common case.
'''

#FIXME: Enforce BETSE initialization on the first import of any BETSE submodule
#by shifting the entirety of this submodule to "betse.__init__" *AFTER*
#stripping out the "betse.lib" and "betse.util" subpackages into a new
#third-party "brutil" package. Why after? Because shifting this submodule now
#would prevent downstream consumers (e.g., BETSEE) from setting the application
#metadata singleton to a downstream preference, which would be bad. Indeed,
#this is the best justification for "brutil" we've stumbled over yet.

# ....................{ IMPORTS                           }....................
from betse.util.app.meta import appmetaone as _appmetaone

# ....................{ MAIN                              }....................
# Instantiate and set a BETSE-specific application metadata singleton if the
# appmetaone.set_app_meta() function has yet to be called elsewhere.
_app_meta = _appmetaone.set_app_meta_betse_if_unset()

# Initialize all mandatory third-party dependencies if the
# _app_meta.init_libs() method has yet to be called elsewhere.
_app_meta.init_libs_if_needed()

# ....................{ CLEANUP                           }....................
# Delete *ALL* attributes (including callables) defined above, preventing the
# package namespace from being polluted with these attributes.
del _appmetaone, _app_meta
