#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level logging configuration functionality.
'''

# ....................{ IMPORTS                           }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To avoid circular import dependencies, avoid importing from *ANY*
# application-specific modules at the top-level -- excluding those explicitly
# known *NOT* to import from this module. Since all application-specific
# modules should *ALWAYS* be able to safely import from this module at any
# scoping level, circularities are best avoided here rather than elsewhere.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

from betse.exceptions import BetseLogException
from betse.util.type.types import type_check

# ....................{ GLOBALS                           }....................
_log_conf = None
'''
Singleton logging configuration for the current Python process.

This configuration provides access to root logger handlers. In particular, this
simplifies modification of logging levels at runtime (e.g., in response to
command-line arguments or configuration file settings).
'''

# ....................{ EXCEPTIONS                        }....................
def die_unless_log_conf() -> None:
    '''
    Raise an exception unless a logging configuration already exists.

    Equivalently, this function raises an exception if *no* logging
    configuration currently exists.

    Raises
    ----------
    BetseLogException
        If no singleton logging configuration currently exists, typically due
        to this function having been called either:

        * *Before* the startup-time call to the :func:`init` function.
        * *After* the shutdown-time call to the :func:`deinit` function.
    '''

    # If no logging configuration currently exists, raise an exception.
    if _log_conf is None:
        raise BetseLogException(
            'Logging unconfigured (e.g., due to either '
            'logconf.init() not having been called or '
            'logconf.deinit() having been called).'
        )

# ....................{ INITIALIZERS                      }....................
@type_check
def init() -> None:
    '''
    Enable our default logging configuration for the active Python process if
    this configuration has yet to be enabled *or* reduce to a noop otherwise
    (e.g., if this method has already been called).

    Specifically, this function instantiates the private :data:`_log_conf`
    singleton to an instance of the application-specific :class:`LogConf`
    class, publicly accessible via the :func:`get_log_conf` module getter. This
    singleton defines sane default filters, formatters, and handlers for the
    root logger, which callers may customize by setting singleton properties.
    '''

    # Avoid circular import dependencies.
    from betse.util.io.log import logs
    from betse.util.io.log.conf.logconfcls import LogConf

    # Module-scoped variables to be set below.
    global _log_conf

    # If a logging configuration already exists...
    if _log_conf is not None:
        # Log a non-fatal warning.
        logs.log_warning(
            'Logging already configured (e.g., due to '
            'logconf.init() having been called).'
        )

        # Reduce to a noop.
        return

    # Instantiate this singleton global with the requisite defaults.
    # print('Reinitializing logging.')
    _log_conf = LogConf()

    # Log this initialization *AFTER* guaranteeing logging sanity.
    logs.log_debug('Initialized singleton logging configuration.')


@type_check
def deinit() -> None:
    '''
    Disable our default logging configuration for the active Python process if
    this configuration has yet to be disabled *or* reduce to a noop otherwise
    (e.g., if this method has already been called).

    Specifically, this function calls the :meth:`LogConf.deinit` method of the
    private :data:`_log_conf` singleton and then nullifies that singleton.
    '''

    # Avoid circular import dependencies.
    from betse.util.io.log import logs

    # Module-scoped variables to be set below.
    global _log_conf

    # If no logging configuration exists, silently noop.
    #
    # Note that this common edge case occurs when an exception is raised early
    # at application startup *BEFORE* the associated init() function is called.
    # Ergo, this constitutes neither a fatal error nor non-fatal warning.
    if _log_conf is None:
        return
    # Else, a logging configuration currently exists.

    # Log this deinitialization *AFTER* guaranteeing logging sanity.
    logs.log_debug('Deinitializing singleton logging configuration...')

    # Deinitialize this logging configuration.
    _log_conf.deinit()

    # Nullify this singleton global for safety *AFTER* all other actions above.
    _log_conf = None

# ....................{ GETTERS                           }....................
def get_log_conf() -> 'LogConf':
    '''
    Singleton logging configuration for the active Python process.

    Raises
    ----------
    BetseLogException
        If no singleton logging configuration currently exists.
    '''

    # If no logging configuration exists, raise an exception.
    die_unless_log_conf()

    # Else, a logging configuration exists. Return this configuration.
    return _log_conf


def get_metadata():
    '''
    Ordered dictionary synopsizing the current logging configuration.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.iterable.mapping.mapcls import OrderedArgsDict

    # Return this dictionary.
    return OrderedArgsDict(
        'file', _log_conf.filename,
        'file level', _log_conf.file_level.name.lower(),
        'verbose', str(_log_conf.is_verbose).lower(),
    )
