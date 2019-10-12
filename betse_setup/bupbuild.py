#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level custom :mod:`setuptools`-specific :class:`ScriptWriter` monkey patch.
'''

# ....................{ IMPORTS                           }....................
from distutils.errors import DistutilsClassError
from setuptools.command.easy_install import ScriptWriter

# ....................{ GLOBALS                           }....................
_SCRIPTWRITER_GET_ARGS_OLD = None
'''
Original (i.e., pre-monkey-patched) implementation of the
:meth:`ScriptWriter.get_args` class method, which our monkey-patch
implementation conditionally calls as needed.
'''

# ....................{ ADDERS                            }....................
def add_subcommand(setup_options: dict, custom_metadata: dict) -> None:

    # Globals to be defined below.
    global _SCRIPTWRITER_GET_ARGS_OLD

    # print(
    #     'Monkey-patching class method '
    #     'setuptools.command.easy_install.ScriptWriter.get_args()...')

    # If this install of setuptools does *NOT* define a "ScriptWriter" class
    # defining the subsequently monkey-patched class method, this install is
    # either broken *OR* of an unsupported version. In either case, raise an
    # exception.
    if not hasattr(ScriptWriter, 'get_args'):
        raise DistutilsClassError(
            'Class method '
            'setuptools.command.easy_install.ScriptWriter.get_args() not '
            'found. The current version of setuptools is either broken '
            '(unlikely) or unsupported (likely).'
        )

    # Preserve the existing implementation of this class method, which our
    # monkey-patch implementation conditionally calls as needed.
    _SCRIPTWRITER_GET_ARGS_OLD = ScriptWriter.get_args

    # Monkey-patch this class method.
    ScriptWriter.get_args = _scriptwriter_get_args_patched

# ....................{ PATCHES                           }....................
# Since this submodule is designed to be copy-and-pasted verbatim into
# downstream consumers (e.g., BETSEE), this implementation is merely a proxy
# that makes no assumptions on whether or not BETSE is currently installed.
@classmethod
def _scriptwriter_get_args_patched(cls, *args, **kwargs):
    '''
    Monkey-patched :meth:`ScriptWriter.get_args` class method proxy.

    This proxy defers to the original implementation of this method until
    :mod:`setuptools` installs BETSE, at which point this proxy defers to
    BETSE's monkey patch of this method.
    '''

    # Attempt to import BETSE.
    try:
        import betse
    # If setuptools has yet to install BETSE, defer to the original
    # implementation of this method.
    except ImportError:
        yield from _SCRIPTWRITER_GET_ARGS_OLD(*args, **kwargs)
        return
    # Else, setuptools has already installed BETSE. In this case, apply our
    # BETSE-specific monkey patch.

    # Defer importation of dependencies that may have yet to be installed.
    from betse import metadata
    from betse.lib.setuptools.command import supcmdbuild

    # Properly monkey-patch the ScriptWriter.get_args() class method.
    supcmdbuild.init(
        package_names={metadata.PACKAGE_NAME,},
        scriptwriter_get_args_old=_SCRIPTWRITER_GET_ARGS_OLD,
    )

    # Validate that this method differs from that method, preventing
    # unwanted infinite recursion.
    if ScriptWriter.get_args is _scriptwriter_get_args_patched:
        raise DistutilsClassError(
            'Class method '
            'setuptools.command.easy_install.ScriptWriter.get_args() not '
            'patched by supcmdbuild.init().'
        )

    # Defer to that method.
    yield from ScriptWriter.get_args(*args, **kwargs)
