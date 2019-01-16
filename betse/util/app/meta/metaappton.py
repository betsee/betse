#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2017-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level **application metadata singleton** (i.e., application-wide object
synopsizing application metadata via read-only properties).
'''

#FIXME: Raise an exception if running with superuser privelages. To avoid
#inevitable security issues, Under *NO* circumstances should either BETSE or
#BETSEE ever be run with elevated permissions. To do so, see the following
#canonical answer:
#    https://stackoverflow.com/a/1026626/2809027

#FIXME: Print memory-related metadata when running "betse info" *AND* log
#non-fatal warnings when BETSE is run under a low-memory environment (e.g.,
#< 4GB available free memory). To do so, note the following canonical API:
#    psutil.Process(os.getpid()).get_memory_info()

#FIXME: The "~/.betse" directory grows fairly large fairly quickly. It'd be
#great to emit non-fatal warnings if its size exceeds some reasonable threshold
#(e.g., 1MB).

# ....................{ IMPORTS                           }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To avoid race conditions during setuptools-based installation, this
# module may import *ONLY* from modules guaranteed to exist at the start of
# installation. This includes all standard Python and application modules but
# *NOT* third-party dependencies, which if currently uninstalled will only be
# installed at some later time in the installation. Likewise, to avoid circular
# import dependencies, the top-level of this module should avoid importing
# application modules where feasible.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

from betse.exceptions import BetseMetaAppException
from betse.util.io.log import logs
from betse.util.type.types import type_check

# ....................{ GLOBALS                           }....................
_app_meta = None
'''
**Application metadata singleton** (i.e., application-wide object synopsizing
application metadata via read-only properties).

Caveats
----------
For safety, callers are advised to call the :func:`get_app_meta` getter safely
returning this private singleton rather than directly accessing this private
singleton unsafely.
'''

# ....................{ EXCEPTIONS                        }....................
def die_if_app_meta() -> None:
    '''
    Raise an exception if the application metadata singleton has been set
    (e.g., by a prior call to the :func:`set_app_meta` function).

    Raises
    ----------
    BetseMetaAppException
        If this singleton has already been set.

    See Also
    ----------
    :func:`is_app_meta`
        Further details.
    '''

    # If an application metadata singleton exists, raise an exception.
    if is_app_meta():
        raise BetseMetaAppException(
            'Application metadata singleton already defined '
            '(i.e., metaappton.set_app_meta() already called).')


def die_unless_app_meta() -> None:
    '''
    Raise an exception unless the application metadata singleton has been set
    (e.g., by a prior call to the :func:`set_app_meta` function).

    Equivalently, this function raises an exception if this singleton has *not*
    yet been set.

    Raises
    ----------
    BetseMetaAppException
        If this singleton has *not* yet been set.

    See Also
    ----------
    :func:`is_app_meta`
        Further details.
    '''

    # If no application metadata singleton exists, raise an exception.
    if not is_app_meta():
        raise BetseMetaAppException(
            'Application metadata singleton undefined '
            '(i.e., metaappton.set_app_meta() not called).')

# ....................{ TESTERS                           }....................
def is_app_meta() -> bool:
    '''
    ``True`` only if the application metadata singleton has been set (e.g., by
    a prior call to the :func:`set_app_meta` function).
    '''

    return _app_meta is not None

# ....................{ GETTERS                           }....................
# Avoid circular import dependencies.
def get_app_meta() -> 'betse.util.app.meta.metaappabc.MetaAppABC':
    '''
    **Application metadata singleton** (i.e., application-wide object
    synopsizing application metadata via read-only properties) if this
    singleton has already been set by a prior call to the :func:`set_app_meta`
    function *or* raise an exception otherwise (i.eg., if that function has yet
    to be called).

    Returns
    ----------
    MetaAppABC
        Application metadata singleton defined by the most recent call to the
        :func:`set_app_meta` function.

    Raises
    ----------
    BetseMetaAppException
        If the :func:`set_app_meta` function has yet to be called.
    '''

    # If no application metadata singleton has been set, raise an exception.
    die_unless_app_meta()
    # Else, an application metadata singleton has been set.

    # Return this sisgleton.
    return _app_meta

# ....................{ SETTERS                           }....................
@type_check
def set_app_meta(
    # Avoid circular import dependencies.
    app_meta: 'betse.util.app.meta.metaappabc.MetaAppABC') -> None:
    '''
    Set the **application metadata singleton** (i.e., application-wide object
    synopsizing application metadata via read-only properties) if this function
    has not already been called *or* raise an exception otherwise (i.e., if
    this function has already been called).

    Caveats
    ----------
    **This function is not intended to be called explicitly.** While callers
    may safely do so, doing so should be entirely redundant. Why? Because the
    :meth:`betse.util.app.meta.metaappabc.MetaAppABC.__init__` method already
    does implicitly at instantiation time.

    **This function intentionally performs no logging.** Doing so would be
    unproductive. The first call to this function is implicitly performed by
    the :func:`betse.util.app.meta.metaappabc.MetaAppABC.__init__` method
    *before* logging has been configured. All logging performed by that call
    (but *not* subsequent calls) would be silently squelched, which any sane
    caller would interpret to be a bug.

    Parameters
    ----------
    app_meta : MetaAppABC
        Application metadata singleton to be set.
    '''

    # Enable this singleton global to be overwritten be the passed parameter.
    global _app_meta

    # If this singleton has already been set, raise an exception.
    die_if_app_meta()

    # Set this singleton global to this caller-specific singleton.
    _app_meta = app_meta

# ....................{ UNSETTERS                         }....................
@type_check
def unset_app_meta() -> None:
    '''
    Unset the **application metadata singleton** (i.e., application-wide object
    synopsizing application metadata via read-only properties).

    Equivalently, this function resets this singleton to its initial state
    (i.e., ``None``).
    '''

    # Enable this singleton global to be overwritten be the passed parameter.
    global _app_meta

    # Log this attempt.
    logs.log_debug('Unsetting application metadata singleton...')

    # Revert this singleton global to its initial state.
    _app_meta = None

# ....................{ MAKERS                            }....................
# Avoid circular import dependencies.
def make_app_meta_betse(*args, **kwargs) -> (
    'betse.util.app.meta.metaappabc.MetaAppABC'):
    '''
    Instantiate and set a BETSE-specific application metadata singleton if the
    :func:`set_app_meta` function has not already been called *and*, in either
    case, (re)initialize this singleton

    This is a convenience function simplifying BETSE initialization for
    low-level edge-case automation (e.g.,
    :mod:`betse.science.__init__`, tests). This function is specific to BETSE
    and hence inappropriate for general-purpose application initialization, as
    required by downstream consumers (e.g., the BETSEE GUI).

    Caveats
    ----------
    **This function does not initialize mandatory third-party dependencies.**
    To permit callers to configure such initialization, callers are required to
    explicitly call the
    :meth:`betse.util.app.meta.metaappabc.MetaAppABC.init_libs` method on the
    object returned by this function.

    Parameters
    ----------
    All parameters are passed as is to the
    :meth:`betse.util.app.meta.metaappabc.MetaAppABC.init_sans_libs` method.

    Returns
    ----------
    betse.util.app.meta.metaappabc.MetaAppABC
        Application metadata singleton.
    '''

    # Avoid circular import dependencies.
    from betse.metaapp import BetseMetaApp

    # If no application metadata singleton has been instantiated, do so. Note
    # that doing so implicitly calls the metaappton.set_app_meta() function on
    # our behalf, which is certainly nice.
    if not is_app_meta():
        BetseMetaApp(*args, **kwargs)
    # An application metadata singleton has now been instantiated.

    # Return this singleton.
    return get_app_meta()
