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

# ....................{ IMPORTS                           }....................
from betse.util.app.meta import metaappton as _metaappton

# ....................{ MAIN                              }....................
# Instantiate and set a BETSE-specific application metadata singleton if the
# metaappton.set_app_meta() function has yet to be called elsewhere.
_app_meta = _metaappton.make_app_meta_betse()

# Initialize all mandatory third-party dependencies if the
# _app_meta.init_libs() method has yet to be called elsewhere.
_app_meta.init_libs_if_needed()

# ....................{ CLEANUP                           }....................
# Delete *ALL* attributes (including callables) defined above, preventing the
# package namespace from being polluted with these attributes.
del _metaappton, _app_meta
