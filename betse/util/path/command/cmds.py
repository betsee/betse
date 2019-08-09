#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level **command** (i.e., external executable file) facilities.
'''

# ....................{ IMPORTS                           }....................
import sys
from betse.exceptions import BetseCommandException
from betse.util.type.decorator.decmemo import func_cached
from betse.util.type.types import type_check, StrOrNoneTypes

# ....................{ EXCEPTIONS                        }....................
@type_check
def die_unless_command(filename: str, reason: StrOrNoneTypes = None) -> None:
    '''
    Raise an exception unless a command with the passed filename exists.

    Parameters
    ----------
    filename : str
        Either the basename *or* the absolute or relative filename of the
        executable file to be validated.
    reason : optional[str]
        Human-readable sentence fragment to be embedded in this exception's
        message (e.g., ``due to "pyside2-tools" not being installed``).
        Defaults to ``None``, in which case this message has no such reason.

    Raises
    ----------
    BetseCommandException
        If this command does *not* exist.

    See Also
    ----------
    :func:`is_command`
        Further details.
    '''

    # If this command does *NOT* exist...
    if not is_command(filename):
        # Exception message to be raised.
        message = 'Command "{}" not found'.format(filename)

        # If an exception reason was passed, embed this reason in this message.
        if reason is not None:
            message += ' ({})'.format(reason)

        # Finalize this message.
        message += '.'

        # Raise this exception.
        raise BetseCommandException(message)

# ....................{ TESTERS                           }....................
@type_check
def is_command(filename: str) -> bool:
    '''
    ``True`` only if a command with the passed filename exists.

    This is the case if this path is either:

    * The basename of an executable file in the current ``${PATH}``.
    * The relative or absolute path of an executable file.

    Parameters
    ----------
    pathname : str
        Either the basename *or* the absolute or relative filename of the
        executable file to be tested.

    Returns
    ----------
    bool
        ``True`` only if this command exists.
    '''

    # Avoid circular import dependencies.
    from betse.util.path import files, pathnames
    from betse.util.path.command import cmdpath

    # This path is that of an existing command if and only if either...
    return (
        # This path is that of an executable file *OR*
        files.is_executable(filename) or (
        # This path is that of a basename in the current ${PATH}.
        pathnames.is_basename(filename) and cmdpath.is_pathable(filename))
    )

# ....................{ GETTERS                           }....................
@func_cached
def get_current_basename() -> str:
    '''
    Basename of the command originating the active Python interpreter.

    If this interpreter is interpreting a block of arbitrary runtime code
    passed to this interpreter on the command line via Python's ``-c`` option
    (e.g., due to being called by a distributed ``pytest-xdist`` test), this
    function unconditionally returns the basename of the current application
    (e.g., ``betse``) rather than ``-c``. Why? We can all agree that ``-c`` is
    an unexpected and rather non-human-readable basename, which is bad.
    '''

    # Avoid circular import dependencies.
    from betse.util.app.meta import appmetaone
    from betse.util.path import pathnames

    # Raw absolute or relative path of the current command.
    current_basename = sys.argv[0]

    # If this is the non-human-readable "-c" Python interpreter option,
    # substitute this with the human-readable basename of this application.
    if current_basename == '-c':
        current_basename = appmetaone.get_app_meta().package_name
    # Else, reduce this absolute or relative path to a basename.
    else:
        current_basename = pathnames.get_basename(current_basename)

    # Return this basename.
    return current_basename
