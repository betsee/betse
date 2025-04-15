#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2025 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Unit tests for the :mod:`betse.util.path.dirs` submodule.
'''

# ....................{ IMPORTS                            }....................
from betse.util.test.pytest.mark.pytskip import (
    # skip_if_ci_gitlab,
    skip_if_os_windows_vanilla,
)

# ....................{ TESTS                              }....................
def test_dirs_get_mtime_newest(
    betse_temp_dir: 'py._path.local.LocalPath') -> None:
    '''
    Unit test the :func:`betse.util.path.dirs.get_mtime_recursive_newest` and
    related :func:`betse.util.path.paths.get_mtime_recursive_newest` functions
    with a directory tree fabricated in the passed temporary directory.

    Parameters
    ----------
    betse_temp_dir : py._path.local.LocalPath
        Object encapsulating a temporary directory isolated to this test.
    '''

    # Defer heavyweight imports.
    from betse.util.path import dirs, paths

    # Absolute paths of a subdirectory and file with arbitrary basenames
    # residing in this temporary directory.
    subdirpath  = betse_temp_dir.join("Quarks_Bar")
    subfilepath = betse_temp_dir.join("Grill")

    # Absolute paths of a subdirectory and two files with arbitrary basenames
    # residing in the above subdirectory residing in this temporary directory.
    subsubdirpath   = subdirpath.join("Gaming_House")
    subsubfilepath1 = subdirpath.join("And")
    subsubfilepath2 = subdirpath.join("Holosuite_Arcade")

    # Create these subdirectories and files, ensuring (...get it, "ensuring"?)
    # that the last such path is that asserted to be the most recent below.
    subdirpath.ensure(dir=True)
    subfilepath.ensure(file=True)
    subsubdirpath.ensure(dir=True)
    subsubfilepath1.ensure(file=True)
    subsubfilepath2.ensure(file=True)

    # Relevant paths as strings rather than "LocalPath" objects.
    dirname = str(betse_temp_dir)
    subsubdirname = str(subsubdirpath)
    subsubfilename2 = str(subsubfilepath2)

    # Ensure the last path created above to be the most recent.
    #
    # Note that *ALL* mtimes (i.e., modification times) here are intentionally
    # truncated from floating-point numbers to integers, implicitly discarding
    # the fractional parts of these numbers. Most modern filesystems maintain
    # fine-grained nanosecond resolution on path timestamps. When a directory
    # tree comprising one parent directory and one or more child subdirectories
    # are recursively created (e.g., via the ensure() method called above), the
    # timestamps associated with these directories often vary with respect to
    # fine-grained nanosecond resolution. Previously, this unit test ignored
    # this distinction and attempted to assert strict equality between these
    # fine-grained timestamps. Doing so raised exceptions resembling:
    #    >       assert (
    #                dirs.get_mtime_recursive_newest(dirname) ==
    #                paths.get_mtime_nonrecursive(subsubfilename2)
    #            )
    #    E       AssertionError: assert 1653971645.6776903 == 1653971645.677683
    #
    # The above output suggests that modern filesystems can reliably provide at
    # most four fractional digits of resolution on path timestamps. For
    # portability, this test goes one step further and simply assumes *ALL*
    # fractional digits of path timestamps to be unreliable. *shakes head, fam*
    assert (
        int(dirs.get_mtime_recursive_newest(dirname)) ==
        int(paths.get_mtime_nonrecursive(subsubfilename2))
    )

    # Update the mtime of an arbitrary subdirectory.
    subsubdirpath.setmtime()

    # Ensure this subdirectory to now be the most recent.
    assert (
        int(dirs.get_mtime_recursive_newest(dirname)) ==
        int(paths.get_mtime_nonrecursive(subsubdirname))
    )

    # Ensure this subdirectory to now be the most recent when queried through
    # the related paths.get_mtime_recursive_newest() function.
    assert (
        int(paths.get_mtime_recursive_newest(
            (dirname, subsubdirname, subsubfilename2))) ==
        int(paths.get_mtime_nonrecursive(subsubdirname))
    )


#FIXME: Resolve the GitLab-CI-specific issue described under "Caveats" below.
#Ideally, this and *ALL* tests should run as is under *ALL* CI hosts. (The
#underlying culprit probably a recent change to GitLab's internal management of
#git repositories - possibly governed by the ${GIT_STRATEGY} variable.)
# @skip_if_ci_gitlab()
@skip_if_os_windows_vanilla()
def test_packages_init() -> None:
    '''
    Unit test the :func:`betse.util.path.dirs.recurse_subdirnames` function by
    validating that all **non-data subdirectories** (i.e., subdirectories
    containing only pure-Python) of *all* top-level package directories of this
    project (i.e., :mod:`betse` and :mod:`betse_test`) contain the mandatory
    ``__init__.py`` special file.

    Caveats
    -------
    **This test is currently incompatible with Microsoft Windows.** Why? Windows
    filesystems are case-insensitive in such a way as to inhibit the sanity of
    this test. Note the differing case in these two local variables of this test
    under a recent GitHub Actions-based Windows test workflow:

        # Note the "Lib" versus "lib" subdirectory basenames in these variables.
        EXCLUDE_DIRNAMES = {
            'D:\\a\\betse\\betse\\.tox\\py310-coverage\\Lib\\site-packages\\betse\\data',
            'D:\\a\\betse\\betse\\betse_test\\_data'
        }
        package_subdirname = (
            'D:\\a\\betse\\betse\\.tox\\py310-coverage\\lib\\site-packages\\betse\\data')

    **This test is currently incompatible with GitLab-CI's host environment.**
    GitLab-CI recently introduced a ``${GIT_STRATEGY}`` environment variable
    governing the tracking of git repositories. Previously, GitLab-CI appears
    to have performed a full-fledged clone of repositories; GitLab-CI appears
    to now be performing incremental pulls of repositories, instead. While
    ostensibly optimal from the efficiency perspective, the latter approach is
    clearly suboptimal from the reproducibility perspective. Each additional
    pull of a repository introduces incremental artifacts into the local
    copy of that repository, whose contents then begin to erroneously diverge.

    Why? Since ``git`` tracks files rather than directories, incremental pulls
    fail to implicitly prune formerly non-empty directories that have since
    been rendered empty by the subsequent removal of all files previously
    contained within those directories.

    While typically innocuous, the existence of empty subdirectories triggers
    false negatives from otherwise working functional and unit tests. In this
    case, this unit test recursively validates that all application
    (sub)packages define a mandatory ``__init__.py`` script, thus failing on
    the first empty non-data subdirectory. Until the underlying host-specific
    issue is resolved, this test is necessarily ignored on this host.
    '''

    # Defer heavyweight imports.
    import betse, betse_test
    from betse.util.app.meta import appmetaone
    from betse.util.io.log import logs
    from betse.util.path import dirs, files, pathnames
    from betse.util.py.module import pymodule
    from betse.util.type.text.string import strs

    # Tuple of all top-level packages.
    PACKAGES = (betse, betse_test)

    # Tuple of the absolute dirnames of all top-level package directories.
    PACKAGE_DIRNAMES = (pymodule.get_dirname(package) for package in PACKAGES)

    # Set of the absolute dirnames of all top-level directories to be .
    EXCLUDE_DIRNAMES = {
        appmetaone.get_app_meta().data_dirname,
        appmetaone.get_app_meta().test_data_dirname,
    }

    # For each such dirname...
    for package_dirname in PACKAGE_DIRNAMES:
        # Log this inspection.
        logs.log_info('')
        logs.log_info(
            'Searching "%s/" for non-package subdirectories...',
            pathnames.get_basename(package_dirname))

        # For the absolute direname of each direct and transitive subdirectory
        # of this package directory...
        for package_subdirname in dirs.recurse_subdirnames(package_dirname):
            # If this is either (in decreasing order of test efficiency)...
            #
            if (
                # A cache subdirectory *OR*...
                pathnames.get_basename(package_subdirname) == '__pycache__' or
                # An empty subdirectory *OR*...
                dirs.is_empty(package_subdirname) or
                # Either an excludable directory or subdirectory of such
                # directory...
                any(
                    strs.is_prefix(
                        text=package_subdirname, prefix=exclude_dirname)
                    for exclude_dirname in EXCLUDE_DIRNAMES
                )
            # ...then this subdirectory is guaranteed to contain no subpackages
            # and hence be safely ignorable.
            ):
                # Log this exclusion.
                logs.log_info('Excluding subdirectory: %s', package_subdirname)

                # Continue to the next subdirectory.
                continue
            # Else, this subdirectory is expected to be a subpackage.

            # Log this inspection.
            logs.log_info('Including subdirectory: %s', package_subdirname)

            # If this subdirectory does *NOT* contain an "__init__.py" file,
            # raise an exception.
            files.join_or_die(package_subdirname, '__init__.py')
