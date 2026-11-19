#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2025 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Unit tests for the :mod:`betse.util.path.files` submodule.
'''

# ....................{ IMPORTS                            }....................
from pytest import raises

# ....................{ TESTS                              }....................
def test_files_die_unless_file_pass(
    betse_temp_dir: 'py._path.local.LocalPath') -> None:
    '''
    Unit test the :func:`betse.util.path.files.die_unless_file` validator
    against existing files, which should *not* raise exceptions.

    Parameters
    ----------
    betse_temp_dir : py._path.local.LocalPath
        Object encapsulating a temporary directory isolated to this test.
    '''

    # Defer heavyweight imports.
    from betse.util.path import files

    # Absolute filenames of two existing files in this temporary directory.
    filepath1 = betse_temp_dir.join('Deep_Space_Nine')
    filepath2 = betse_temp_dir.join('Terok_Nor')
    filepath1.ensure(file=True)
    filepath2.ensure(file=True)

    # Assert this validator to accept these files.
    files.die_unless_file(str(filepath1))
    files.die_unless_file(str(filepath1), str(filepath2))


def test_files_die_unless_file_fail(
    betse_temp_dir: 'py._path.local.LocalPath') -> None:
    '''
    Unit test the :func:`betse.util.path.files.die_unless_file` validator
    against paths that are *not* existing files, which should raise exceptions
    naming the offending path and explaining why that path is invalid.

    Parameters
    ----------
    betse_temp_dir : py._path.local.LocalPath
        Object encapsulating a temporary directory isolated to this test.
    '''

    # Defer heavyweight imports.
    from betse.exceptions import BetseFileException
    from betse.util.path import files

    # Absolute filename of an existing file in this temporary directory.
    filepath_exists = betse_temp_dir.join('Bajor')
    filepath_exists.ensure(file=True)
    filename_exists = str(filepath_exists)

    # Absolute filename of a missing file in this existing directory.
    filename_missing = str(betse_temp_dir.join('Gamma_Quadrant'))

    # Assert this validator to reject this missing file with an exception
    # message naming this file.
    with raises(BetseFileException) as exception_info:
        files.die_unless_file(filename_missing)
    exception_message = str(exception_info.value)
    assert filename_missing in exception_message
    assert 'no such file exists' in exception_message

    # Assert this validator to reject this missing file even when preceded by
    # an existing file, naming the missing rather than the existing file.
    with raises(BetseFileException) as exception_info:
        files.die_unless_file(filename_exists, filename_missing)
    exception_message = str(exception_info.value)
    assert filename_missing in exception_message
    assert filename_exists not in exception_message

    # Absolute dirname of a missing directory in this temporary directory and
    # absolute filename of a missing file in a missing subdirectory of that.
    dirpath_missing = betse_temp_dir.join('Wormhole')
    filename_missing_nested = str(
        dirpath_missing.join('Celestial_Temple').join('Prophets.yaml'))

    # Assert this validator to reject this missing file with an exception
    # message naming the highest missing directory containing this file.
    with raises(BetseFileException) as exception_info:
        files.die_unless_file(filename_missing_nested)
    exception_message = str(exception_info.value)
    assert filename_missing_nested in exception_message
    assert 'directory "{}" does not exist'.format(
        str(dirpath_missing)) in exception_message

    # Assert this validator to reject an existing directory with an exception
    # message explaining that this path is a directory.
    with raises(BetseFileException) as exception_info:
        files.die_unless_file(str(betse_temp_dir))
    assert 'is a directory' in str(exception_info.value)
