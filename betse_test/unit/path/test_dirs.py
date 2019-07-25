#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Unit tests for the :mod:`betse.util.path.dirs` submodule.
'''

# ....................{ IMPORTS                           }....................
from betse.util.test.pytest.mark.pytskip import skip_if_ci_gitlab

# ....................{ TESTS                             }....................
def test_dirs_get_mtime_newest(betse_temp_dir: 'LocalPath') -> None:
    '''
    Unit test the :func:`betse.util.path.dirs.get_mtime_recursive_newest` and
    related :func:`betse.util.path.paths.get_mtime_recursive_newest` functions
    with a directory tree fabricated in the passed temporary directory.

    Parameters
    ----------
    betse_temp_dir : LocalPath
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
    # that the last such path is that asserted to be the most recent below
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
    assert (
        dirs.get_mtime_recursive_newest(dirname) ==
        paths.get_mtime_nonrecursive(subsubfilename2)
    )

    # Update the mtime of an arbitrary subdirectory.
    subsubdirpath.setmtime()

    # Ensure this subdirectory to now be the most recent.
    assert (
        dirs.get_mtime_recursive_newest(dirname) ==
        paths.get_mtime_nonrecursive(subsubdirname)
    )

    # Ensure this subdirectory to now be the most recent when queried through
    # the related paths.get_mtime_recursive_newest() function.
    assert (
        paths.get_mtime_recursive_newest(
            (dirname, subsubdirname, subsubfilename2)) ==
        paths.get_mtime_nonrecursive(subsubdirname)
    )


#FIXME: Resolve the GitLab-CI-specific issue described under "Caveats" below.
#Ideally, this and *ALL* tests should run as is under *ALL* CI hosts. (The
#underlying culprit probably a recent change to GitLab's internal management of
#git repositories - possibly governed by the ${GIT_STRATEGY} variable.)
@skip_if_ci_gitlab()
def test_packages_init() -> None:
    '''
    Unit test the :func:`betse.util.path.dirs.recurse_subdirnames` function by
    validating that all **non-data subdirectories** (i.e., subdirectories
    containing only pure-Python) of *all* top-level package directories of this
    project (i.e., :mod:`betse`, :mod:`betse_setup`, :mod:`betse_test`) contain
    the mandatory ``__init__.py`` special file.

    Caveats
    ----------
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
    import betse, betse_setup, betse_test
    from betse.util.app.meta import appmetaone
    from betse.util.io.log import logs
    from betse.util.path import dirs, files, pathnames
    from betse.util.py.module import pymodule
    from betse.util.type.text.string import strs

    # Tuple of all top-level packages.
    PACKAGES = (betse, betse_setup, betse_test)

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
