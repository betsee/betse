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
from betse.metaapp import BetseMetaApp
from betse.util.app.meta import metaappton

# ....................{ INITIALIZERS                      }....................
# If no application metadata singleton has been instantiated, do so. Note that
# doing so already calls the metaappton.set_app_meta() function, which is nice.
if not metaappton.is_app_meta():
    BetseMetaApp()
# An application metadata singleton has now been instantiated.

# Initialize all mandatory dependencies with sane defaults. For safety, do so
# regardless of whether an application metadata singleton has been
# instantiated. Since the MetaAppABC.__init__() method does *NOT* implicitly
# call the MetaAppABC.init_libs() method, do so explicitly here to guarantee
# this application to be fully initialized.
metaappton.get_app_meta().init_libs()

# ....................{ CLEANUP                           }....................
# Delete *ALL* attributes (including callables) defined above, preventing the
# package namespace from being polluted with these attributes.
del BetseMetaApp, metaappton
