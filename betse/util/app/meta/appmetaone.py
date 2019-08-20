#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
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
    Raise an exception if the application metadata singleton already exists
    (e.g., due to a prior call to the :func:`set_app_meta` function).

    Raises
    ----------
    BetseMetaAppException
        If this singleton already exists.

    See Also
    ----------
    :func:`is_app_meta`
        Further details.
    '''

    # If an application metadata singleton exists, raise an exception.
    if is_app_meta():
        raise BetseMetaAppException(
            'Application metadata singleton already defined '
            '(i.e., appmetaone.set_app_meta() already called).')


def die_unless_app_meta() -> None:
    '''
    Raise an exception unless an application metadata singleton exists (e.g.,
    due to a prior call to the :func:`set_app_meta` function).

    Equivalently, this function raises an exception if this singleton does
    *not* exist.

    Raises
    ----------
    BetseMetaAppException
        If this singleton does *not* exist.

    See Also
    ----------
    :func:`is_app_meta`
        Further details.
    '''

    # If no application metadata singleton exists, raise an exception.
    if not is_app_meta():
        raise BetseMetaAppException(
            'Application metadata singleton undefined '
            '(i.e., appmetaone.set_app_meta() not called).')

# ....................{ TESTERS                           }....................
def is_app_meta() -> bool:
    '''
    ``True`` only if an application metadata singleton exists (e.g., from a
    prior call to the :func:`set_app_meta` function).
    '''

    return _app_meta is not None

# ....................{ GETTERS                           }....................
# Avoid circular import dependencies.
def get_app_meta() -> 'betse.util.app.meta.appmetaabc.AppMetaABC':
    '''
    **Application metadata singleton** (i.e., application-wide object
    synopsizing application metadata via read-only properties) if this
    singleton has already been set by a prior call to the :func:`set_app_meta`
    function *or* raise an exception otherwise (i.eg., if that function has yet
    to be called).

    Returns
    ----------
    AppMetaABC
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
    app_meta: 'betse.util.app.meta.appmetaabc.AppMetaABC') -> None:
    '''
    Set the **application metadata singleton** (i.e., application-wide object
    synopsizing application metadata via read-only properties) if this is the
    first call to this function *or* raise an exception otherwise (i.e., if
    this function has already been called).

    Caveats
    ----------
    **This function intentionally performs no logging.** Doing so would be
    unproductive. The first call to this function is implicitly performed by
    the :func:`betse.util.app.meta.appmetaabc.AppMetaABC.__init__` method
    *before* logging has been configured. All logging performed by that call
    (but *not* subsequent calls) would be silently squelched, which any sane
    caller would interpret to be a bug.

    Parameters
    ----------
    app_meta : AppMetaABC
        Application metadata singleton to be set.

    Raises
    ----------
    BetseMetaAppException
        If this singleton has already been set.

    See Also
    ----------
    :meth:`betse.util.app.meta.appmetaabc.AppMetaABC.__init__`
        Higher-level method encapsulating this lower-level function.
    '''

    # Enable this singleton global to be overwritten be the passed parameter.
    global _app_meta

    # If this singleton has already been set, raise an exception.
    die_if_app_meta()

    # Set this singleton global to this caller-specific singleton.
    _app_meta = app_meta

# ....................{ SETTERS ~ unset                   }....................
#FIXME: Preserve but refactor as follows:
#
#* Define a new set_app_meta_if_unset() function in this submodule with the
#  following signature, where "cls" is the "AppMetaABC" subclass to be
#  conditionally instantiated by this function if needed:
#    def set_app_meta_if_unset(cls: ClassType, *args, **kwargs) -> (
#       'betse.util.app.meta.appmetaabc.AppMetaABC')
#* Refactor set_app_meta_betse_if_unset() in terms of set_app_meta_if_unset().
def set_app_meta_betse_if_unset(*args, **kwargs) -> (
    # Avoid circular import dependencies.
    'betse.util.app.meta.appmetaabc.AppMetaABC'):
    '''
    Set the application metadata singleton to a BETSE-specific singleton (i.e.,
    instance of the :class:`betse.appmeta.BetseAppMeta` class) if the
    :func:`set_app_meta` function has yet to be called *or* reduce to a noop
    otherwise (i.e., if that function has already been called).

    This is a convenience function simplifying BETSE initialization for
    low-level edge-case automation (e.g.,
    :mod:`betse.science.__init__`, tests). This function is specific to BETSE
    and hence inappropriate for general-purpose application initialization, as
    required by downstream consumers (e.g., the BETSEE GUI).

    Caveats
    ----------
    **This function does not initialize mandatory third-party dependencies.**
    To enable callers to configure such initialization, callers are required to
    explicitly call the
    :meth:`betse.util.app.meta.appmetaabc.AppMetaABC.init_libs` method on the
    object returned by this function.

    Parameters
    ----------
    All parameters are passed as is to the
    :meth:`betse.util.app.meta.appmetaabc.AppMetaABC.init_sans_libs` method.

    Returns
    ----------
    betse.util.app.meta.appmetaabc.AppMetaABC
        Application metadata singleton.
    '''

    # Avoid circular import dependencies.
    from betse.appmeta import BetseAppMeta

    # If no application metadata singleton has been instantiated, do so. Note
    # that doing so implicitly calls the appmetaone.set_app_meta() function on
    # our behalf, which is certainly nice.
    if not is_app_meta():
        BetseAppMeta(*args, **kwargs)
    # An application metadata singleton has now been instantiated.

    # Return this singleton.
    return get_app_meta()

# ....................{ UNSETTERS                         }....................
def unset_app_meta() -> None:
    '''
    Unset the **application metadata singleton** (i.e., application-wide object
    synopsizing application metadata via read-only properties) if such a
    singleton exists *or* raise an exception otherwise (i.e., if no such
    singleton exists).

    Equivalently, this function resets this singleton to its initial state
    (i.e., ``None``).

    Raises
    ----------
    BetseMetaAppException
        If this singleton has not yet been set.

    See Also
    ----------
    :func:`deinit`
        Higher-level function encapsulating this lower-level function.
    '''

    # Enable this singleton global to be overwritten be the passed parameter.
    global _app_meta

    # If this singleton has not yet been set, raise an exception.
    die_unless_app_meta()

    # Revert this singleton global to its initial state.
    #
    # Note that this reversion is intentionally *NOT* logged. This method is
    # typically only called as the last operation at application shutdown, at
    # which time the logging configuration has already been deinitialized.
    # While feasible, logging here would do so under Python's default logging
    # configuration in a non-human-readable format upsetting end users: e.g.,
    #     DEBUG:betsee:Unsetting application metadata singleton...
    _app_meta = None

# ....................{ DEINITIALIZERS                    }....................
def deinit() -> None:
    '''
    Deinitialize the **application metadata singleton** (i.e., application-wide
    object synopsizing application metadata via read-only properties) if such a
    singleton exists *or* silently reduce to a noop otherwise (i.e., if no such
    singleton exists).

    This function effectively deinitializes the current application, which
    requires this singleton for rudimentary functionality throughout the
    codebase.

    Caveats
    ----------
    **No application logic may be safely performed after calling this method.**
    This method nullifies this singleton *and* closes open logfile handles,
    both of which are required for basic application logic.
    '''

    # If an application metadata singleton exists, deinitialize this singleton.
    if is_app_meta():
        _app_meta.deinit()
    # Else, no such singleton exists. In this case, *NO* logic (including
    # logging) may be safely performed. Ergo, noop.
