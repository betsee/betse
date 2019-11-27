#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level **process tester** (i.e., utility functions testing and validating
currently running platform-specific processes) facilities.
'''

# ....................{ IMPORTS                           }....................
from betse.exceptions import (
    BetseFunctionUnimplementedException,
    BetseProcessNotFoundException,
)
from betse.util.type.types import type_check

# ....................{ TESTERS ~ parent                  }....................
#FIXME: This tester is currently uncalled. Before calling this tester, the
#optional "psutil" dependency *MUST* be uncommented in "betse.metadeps".
@type_check
def is_parent_command(command_basename: str) -> bool:
    '''
    ``True`` only if the parent process of the active Python interpreter is
    running an external command with the passed basename optionally suffixed by
    a platform-specific filetype (e.g., ``.exe`` under Windows).

    Parameters
    ----------
    command_basename : str
        Basename of the command to test this process' parent process against.

    Raises
    ------
    BetseOSException
        If the optional :mod:`psutil` dependency is unimportable *and* the
        current platform is neither Linux, macOS, or Windows.
    BetseFunctionUnimplementedException
        If the optional :mod:`psutil` dependency is unimportable.
    BetseProcessNotFoundException
        If either:

        * The optional :mod:`psutil` dependency is importable *and* either:

          * The current process has *no* parent process.
          * The current process had a parent process that has since died.

    Returns
    ----------
    bool
        ``True`` only if the parent process of the active Python interpreter is
        running an external command with the passed basename.
    '''

    # Avoid circular import dependencies.
    from betse.lib import libs
    from betse.util.os import oses
    from betse.util.os.brand import linux, macos, windows

    # If the optional "psutil" dependency is importable, this dependency
    # typically implements more efficient *AND* portable solutions than all
    # alternatives and is hence preferred. In this case, defer to "psutil".
    if libs.is_runtime_optional('psutil'):
        return _is_parent_command_psutil(command_basename)
    # Else if the current platform is Linux, defer to a Linux-specific tester.
    elif linux.is_linux():
        return _is_parent_command_linux(command_basename)
    # Else if the current platform is macOS, defer to a macOS-specific tester.
    elif macos.is_macos():
        return _is_parent_command_macos(command_basename)
    # Else if the current platform is Windows, defer to a Windows-specific
    # tester.
    elif windows.is_windows():
        return _is_parent_command_windows(command_basename)

    # Else, raise an exception.
    oses.die_if_unsupported()


@type_check
def _is_parent_command_psutil(command_basename: str) -> bool:
    '''
    ``True`` only if the parent process of the active Python interpreter is
    running an external command with the passed basename, implemented in terms
    of the optional :mod:`psutil` dependency.

    See Also
    ----------
    https://stackoverflow.com/a/2241047/2809027
        StackOverflow answer strongly inspiring this implementation.
    :func:`is_parent_command`
        Further details.
    '''

    # Avoid circular import dependencies.
    from betse.lib import libs
    from betse.util.os.brand import windows
    from betse.util.path import pathnames

    # Optional "psutil" dependency.
    psutil = libs.import_runtime_optional('psutil')

    # Object encapsulating the current process.
    current_proc = psutil.Process()

    # Object encapsulating the parent process of this process if any *OR*
    # "None" otherwise.
    parent_proc = current_proc.parent()

    # If the current process has *NO* parent process, raise an exception.
    #
    # In theory, all processes except the initial "init" process should
    # *ALWAYS* have a parent process on POSIX-compatible platforms. In
    # practice, "psutil" documentation explicitly states that the
    # psutil.Process.parent() method returns "None" in edge cases. *sigh*
    if parent_proc is None:
        raise BetseProcessNotFoundException(
            'Current process {} parent not found.'.format(current_proc))
    # Else, the current process has a parent process.

    # If...
    if (
        # The current platform is Windows *AND*...
        windows.is_windows() and
        # This command basename is *NOT* suffixed by a filetype...
        not pathnames.is_filetype(command_basename)
    # Then suffix this basename by ".exe". This is required, as Windows
    # executables are *ALWAYS* suffixed by such a filetype.
    ):
        command_basename += '.exe'

    # Attempt to...
    try:
        # Dictionary of all metadata required to implement this tester.
        #
        # Note that this call:
        #
        # * Efficiently retrieves this metadata with one method call rather
        #   than inefficiently distributing this retrieval across multiple
        #   method calls.
        # * Transparently sets the values of all metadata whose retrieval
        #   raises an "AccessDenied" or "ZombieProcess" exception to "None."
        #   Sadly, the "NoSuchProcess" exception is *NOT* handled similarly and
        #   *MUST* thus be explicitly caught below.
        parent_proc_metadata = parent_proc.as_dict(
            attrs=('name', 'exe', 'cmdline'))

        # List of one or more strings comprising the command line invoking
        # this parent process *OR* "None" if retrieving this metadata
        # raised an exception above.
        parent_proc_cmdline = parent_proc_metadata['cmdline']

        # If either...
        return (
            # This command basename is this parent process' name *OR*...
            #
            # Note that this metadata typically requires the least privelages
            # on most platforms *AND* is typically the most efficient such
            # metadata to access. So, this metadata is tested first.
            command_basename == parent_proc_metadata['name'] or
            # This command basename is that of this parent process *OR*...
            command_basename == pathnames.get_basename(
                parent_proc_metadata['exe']) or
            (
                # This parent process' command line was retrievable *AND*...
                parent_proc_cmdline and
                # This command basename is that of the first item of this line.
                command_basename == pathnames.get_basename(
                    parent_proc_cmdline[0])
            )
        )
    # If this parent process died during metadata access, raise an exception.
    except psutil.NoSuchProcess as exception:
        raise BetseProcessNotFoundException(
            'Current process {} parent {} killed.'.format(
                current_proc, parent_proc)
        ) from exception


#FIXME: Implement us up given the StackOverflow answer listed below.
@type_check
def _is_parent_command_linux(command_basename: str) -> bool:
    '''
    ``True`` only if the parent process of the active Python interpreter is
    running an external command with the passed basename, implemented in terms
    of the Linux-specific **process pseudo-filesystem** (i.e., ``/proc``).

    See Also
    ----------
    https://stackoverflow.com/a/24114907/2809027
        StackOverflow answer strongly inspiring this implementation.
    :func:`is_parent_command`
        Further details.
    '''

    raise BetseFunctionUnimplementedException()


#FIXME: Research whether implementing this tester is even feasible.
@type_check
def _is_parent_command_macos(command_basename: str) -> bool:
    '''
    ``True`` only if the parent process of the active Python interpreter is
    running an external command with the passed basename, implemented in terms
    of the macOS-specific... well, who actually knows?

    See Also
    ----------
    :func:`is_parent_command`
        Further details.
    '''

    raise BetseFunctionUnimplementedException()


#FIXME: Incorporate the following brilliant pure-Python solution for parsing
#Windows process metadata without forking external processes:
#    https://stackoverflow.com/a/7110486/2809027
@type_check
def _is_parent_command_windows(command_basename: str) -> bool:
    '''
    ``True`` only if the parent process of the active Python interpreter is
    running an external command with the passed basename, implemented in terms
    of the Windows-specific Win32 API.

    See Also
    ----------
    https://stackoverflow.com/a/7110486/2809027
        StackOverflow answer strongly inspiring this implementation.
    :func:`is_parent_command`
        Further details.
    '''

    raise BetseFunctionUnimplementedException()
