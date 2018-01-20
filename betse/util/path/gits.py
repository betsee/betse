#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level **Git** (i.e., third-party version control system remotely tracking
all changes to this application) facilities.
'''

# ....................{ IMPORTS                            }....................
from betse.exceptions import BetseGitException
from betse.util.io.log.logenum import LogLevel
from betse.util.type.types import type_check

# ....................{ EXCEPTIONS                         }....................
@type_check
def die_unless_worktree(dirname: str) -> None:
    '''
    Raise an exception unless the directory with the passed pathname is a **Git
    working tree** (i.e., contains a ``.git`` subdirectory).

    Parameters
    ----------
    dirname : str
        Relative or absolute pathname of the directory to be validated.

    Raises
    ----------
    BetseGitException
        If this directory is *not* a Git working tree.

    See Also
    ----------
    :func:`is_repo`
        Further details.
    '''

    # If this directory is *NOT* a Git working tree, raise an exception.
    if not is_worktree(dirname):
        raise BetseGitException(
            'Directory "{}" not a Git working tree '
            '(i.e., contains no ".git" subdirectory).'.format(dirname))

# ....................{ TESTERS                            }....................
@type_check
def is_worktree(dirname: str) -> bool:
    '''
    `True` only if the directory with the passed pathname is a **Git
    working tree** (i.e., contains a ``.git`` subdirectory).

    Parameters
    ----------
    dirname : str
        Relative or absolute pathname of the directory to be tested.

    Returns
    ----------
    bool
        `True` only if this directory is a Git working tree.
    '''

    # Avoid circular import dependencies.
    from betse.util.path import dirs, pathnames

    # Absolute or relative pathname of this directory's ".git" subdirectory.
    git_subdirname = pathnames.join(dirname, '.git')

    # Return True only if this subdirectory exists.
    return dirs.is_dir(git_subdirname)

# ....................{ IO                                 }....................
@type_check
def clone_worktree_shallow(
    branch_or_tag_name: str,
    src_dirname: str,
    trg_dirname: str,
) -> None:
    '''
    **Shallowly clone** (i.e., recursively copy without commit history) from the
    the passed branch or tag name of the source Git working tree with the passed
    dirname into the target non-existing directory.

    For space efficiency, the target directory will *not* be a Git working tree
    (i.e., will *not* contain a ``.git`` subdirectory). This directory will only
    contain the contents of the source Git working tree isolated at either:

    * The ``HEAD`` of the branch with the passed name.
    * The commit referenced by the tag with the passed name.

    Caveats
    ----------
    Due to long-standing issues in Git itself, the passed branch or tag name
    *cannot* be the general-purpose SHA-1 hash of an arbitrary commit. While
    performing a shallow clone of an arbitrary commit from such a hash is
    technically feasible, doing so requires running a pipeline of multiple
    sequential Git subprocesses and hence exceeds the mandate of this function.

    Parameters
    ----------
    branch_or_tag_name : str
        Name of the branch or tag to clone from this source Git working tree.
    src_dirname : str
        Relative or absolute pathname of the source Git working tree to clone
        from.
    trg_dirname : str
        Relative or absolute pathname of the target non-existing directory to
        clone into.

    Raises
    ----------
    BetseGitException
        If this source directory is *not* a Git working tree.
    BetseDirException
        If this target directory already exists.

    See Also
    ----------
    https://stackoverflow.com/questions/26135216/why-isnt-there-a-git-clone-specific-commit-option
        StackOverflow question entitled "Why Isn't There A Git Clone Specific
        Commit Option?", detailing Git's current omission of such functionality.
    '''

    # Avoid circular import dependencies.
    from betse.util.path import paths, pathnames
    from betse.util.path.command import cmdrun

    # If this source directory is *NOT* a Git working tree, raise an exception.
    die_unless_worktree(src_dirname)

    # If this target directory already exists, raise an exception.
    paths.die_if_path(trg_dirname)

    # Absolute pathname of the source Git working tree to clone from. The
    # "git clone" command requires "file:///"-prefixed absolute pathnames.
    src_dirname = pathnames.canonicalize(src_dirname)

    # Git-specific URI of this source Git working tree.
    src_dir_uri = 'file://' + src_dirname

    # Tuple of shell words comprising the Git command performing this clone.
    git_command = (
        'git', 'clone',

        # Branch or tag to be cloned.
        '--branch', branch_or_tag_name,

        # Perform a shallow clone.
        '--depth', '1',

        # From this source to target directory.
        src_dir_uri, trg_dirname,
    )

    # Shallowly clone this source to target directory. Contrary to expectation,
    # "git" redirects non-error or -warning output resembling the following to
    # stderr rather than stdout:
    #
    #    Cloning into '/tmp/pytest-of-leycec/pytest-26/cli_sim_backward_compatibility0/betse_old'...
    #    Note: checking out 'd7d6bf6d61ff2b467f9983bc6395a8ba9d0f234e'.
    #
    #    You are in 'detached HEAD' state. You can look around, make experimental
    #    changes and commit them, and you can discard any commits you make in this
    #    state without impacting any branches by performing another checkout.
    #
    #    If you want to create a new branch to retain commits you create, you may
    #    do so (now or later) by using -b with the checkout command again. Example:
    #
    #      git checkout -b <new-branch-name>
    #
    # This is distasteful. Clearly, this output should *NOT* be logged as either
    # an error or warning. Instead, both stdout and stderr output would ideally
    # be logged as informational messages. Unfortunately, the following call
    # appears to erroneously squelch rather than redirect most stderr output:
    #
    #    cmdrun.log_output_or_die(
    #        command_words=git_command,
    #        stdout_log_level=LogLevel.INFO,
    #        stderr_log_level=LogLevel.ERROR,
    #    )
    #
    # Frankly, we have no idea what is happening here -- and it doesn't
    # particularly matter. The current approach, while non-ideal, suffices.
    cmdrun.run_or_die(command_words=git_command)
