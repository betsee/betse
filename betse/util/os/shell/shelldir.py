#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level shell-specific directory facilities.
'''

# ....................{ IMPORTS                            }....................
import os
from betse.util.io.log import logs
from betse.util.type.types import type_check
from contextlib import contextmanager
from types import GeneratorType

# ....................{ GETTERS                            }....................
def get_cwd_dirname() -> str:
    '''
    **Current working dirname** (i.e., absolute path of the current working
    directory (CWD)) of the active Python process.

    Unless subsequently changed, this is the absolute path of the directory from
    which this application was initially run.
    '''

    return os.getcwd()

# ....................{ SETTERS                            }....................
@type_check
def set_cwd(dirname: str) -> None:
    '''
    Set the **current working directory** (CWD) of the active Python process to
    the passed directory.

    This function permanently changes the CWD for the remainder of this process.
    For a robust alternative changing the CWD for a single code block, consider
    using the :func:`setting_cwd` context manager instead.

    Parameters
    -----------
    dirname : str
        Absolute or relative path of the directory to change to.
    '''

    # Log this change.
    logs.log_debug('Changing current working directory to: %s', dirname)

    # Change to this directory.
    os.chdir(dirname)

# ....................{ CONTEXTS                           }....................
@contextmanager
@type_check
def setting_cwd(dirname: str) -> GeneratorType:
    '''
    Context manager setting the **current working directory** (CWD) of the
    active Python process to the passed directory for the duration of this
    context.

    This context manager guaranteeably reverts the CWD to the prior CWD even
    when fatal exceptions are raised (e.g., due to this directory not existing).

    Parameters
    -----------
    dirname : str
        Absolute or relative path of the directory to change to.

    Returns
    -----------
    contextlib._GeneratorContextManager
        Context manager changing the CWD as described above.

    Yields
    -----------
    None
        Since this context manager yields no values, the caller's ``with``
        statement must be suffixed by *no* ``as`` clause.

    See Also
    -----------
    https://stackoverflow.com/a/24176022/2809027
        StackOverflow answer strongly inspiring this implementation.

    Examples
    -----------
    >>> from betse.util.paths import dirs
    >>> print('CWD: ' + dirs.get_cwd_dirname())
    CWD: /home/azrael
    >>> with dirs.current('/home/uriel/urial/nuriel/uryan/jeremiel'):
    ...     print('CWD: ' + dirs.get_cwd_dirname())
    ...     raise ValueError(
    ...         'But unknown, abstracted, brooding secret the dark power hid')
    CWD: /home/uriel/urial/nuriel/uryan/jeremiel
    ValueError: But unknown, abstracted, brooding secret the dark power hid.
    >>> print('CWD: ' + dirs.get_cwd_dirname())
    CWD: /home/azrael
    '''

    # Absolute path of the current CWD.
    dirname_prior = get_cwd_dirname()

    # Temporarily change to the passed directory. Since Python performs this
    # change only if this call raises no exceptions, this call need *NOT* be
    # embedded in the "try" block below.
    set_cwd(dirname)

    # Yield control to the body of the caller's "with" block.
    try:
        yield
    # Revert to the prior CWD even if that block raised an exception.
    finally:
        os.chdir(dirname_prior)
