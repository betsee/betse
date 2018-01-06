#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Unit tests for the :mod:`betse.util.path.dirs` submodule.
'''

# ....................{ IMPORTS                            }....................

# ....................{ TESTS                              }....................
def test_dirs_get_mtime_newest(betse_temp_dir: 'LocalPath') -> None:
    '''
    Unit test the :func:`betse.util.path.dirs.get_mtime_recursive_newest` and
    related :func:`betse.util.path.paths.get_mtime_recursive_newest` functions with a
    directory tree fabricated in the passed temporary directory.

    Parameters
    ----------
    betse_temp_dir : LocalPath
        Object encapsulating a temporary directory isolated to the current test.
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
