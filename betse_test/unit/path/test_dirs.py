#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Unit tests for the :mod:`betse.util.path.dirs` submodule.
'''

# ....................{ IMPORTS                           }....................

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


def test_dirs_recurse_subdirnames() -> None:
    '''
    Unit test the :func:`betse.util.path.dirs.recurse_subdirnames` function by
    validating that all **non-data subdirectories** (i.e., subdirectories
    containing only pure-Python) of *all* top-level package directories of this
    project (i.e., :mod:`betse`, :mod:`betse_setup`, :mod:`betse_test`) contain
    the mandatory ``__init__.py`` special file.
    '''

    # Defer heavyweight imports.
    import betse, betse_setup, betse_test
    from betse import metaapp
    from betse.util.io.log import logs
    from betse.util.py import pymodule
    from betse.util.path import dirs, files, pathnames
    from betse.util.type.text import strs

    # Tuple of all top-level packages.
    PACKAGES = (betse, betse_setup, betse_test)

    # Tuple of the absolute dirnames of all top-level package directories.
    PACKAGE_DIRNAMES = (pymodule.get_dirname(package) for package in PACKAGES)

    # Absolute dirname of this application's top-level data directory.
    DATA_DIRNAME = metaapp.get_app_meta().data_dirname

    # For each such dirname...
    for package_dirname in PACKAGE_DIRNAMES:
        # Log this inspection.
        logs.log_info(
            'Searching "%s/" for non-package subdirectories...',
            pathnames.get_basename(package_dirname))

        # For the absolute direname of each direct and transitive subdirectory
        # of this package directory...
        for package_subdirname in dirs.recurse_subdirnames(package_dirname):
            # If this is either (in decreasing order of efficiency):
            #
            # * A cache subdirectory.
            # * The data directory.
            # * A subdirectory of the data directory.
            # * An empty subdirectory.
            #
            # Then this subdirectory is guaranteed to contain no subpackages
            # and hence be safely ignorable.
            if (
                pathnames.get_basename(package_subdirname) == '__pycache__' or
                strs.is_prefix(text=package_subdirname, prefix=DATA_DIRNAME) or
                dirs.is_empty(package_subdirname)
            ):
                # Log this exclusion.
                logs.log_info('Excluding subdirectory: %s', package_subdirname)

                # Continue to the next subdirectory.
                continue
            # Else, this subdirectory is expected to be a subpackage.

            # Log this inspection.
            logs.log_info('Including subdirectory: %s', package_subdirname)

            # Absolute filename of the "__init__.py" file of this subdirectory.
            package_init_filename = pathnames.join(
                package_subdirname, '__init__.py')

            # If this file does *NOT* exist, raise an exception.
            files.die_unless_file(package_init_filename)
