#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Unit tests for the :mod:`betse.util.path.archives` submodule.
'''

# ....................{ IMPORTS                            }....................
import pytest
from betse_test.util.mark.skip import skip_unless_module

# ....................{ CONSTANTS                          }....................
# Ideally, the existing
# "betse.util.path.archives._ARCHIVE_FILETYPE_TO_MODULE_NAME" dictionary would
# be leveraged to populate this parametrization. Unfortunately, doing so would
# necessitate unsafely importing the "betse.util.path.archives" submodule above,
# ensuring failure of the entire test suite on the inability to import this
# submodule -- absurd and frankly unacceptable test-time fragility. Instead, the
# set of all archive filetypes is redefined in test-oriented manner ignoring
# filetypes whose corresponding optional stdlib modules are unavailable.
#
# See the "_ARCHIVE_FILETYPE_TO_MODULE_NAME" dictionary for further details.
ARCHIVE_FILETYPES = (
    # The higher-level skip_unless_lib_runtime_optional() decorator is
    # intentionally *NOT* called here. Although the following modules could be
    # argued to be optional runtime dependencies of this application requiring
    # addition to the "betse.metadata.RUNTIME_OPTIONAL" dictionary,
    # these modules are pre-packaged with Python itself rather than installed
    # via setuptools. Since these modules declare no "__version__" attribute,
    # there exists no means of validating the satisfiability of these modules as
    # with customary setuptools-installed optional dependencies.
    #
    # The lower-level skip_unless_module() decorator is thus applicable here.
    skip_unless_module('bz2')(('bz2',)),
    skip_unless_module('gzip')(('gz',)),
    skip_unless_module('lzma')(('xz',)),
)
'''
Tuple of all supported archive filetypes, ignoring those whose corresponding
optional stdlib modules are unimportable under the active Python interpreter and
hence *not* installed at Python installation time.
'''


ARCHIVE_BYTES = bytes(
    (
        "The economics of the future is somewhat different. "
        "You see, money doesn't exist in the 24th century. "
        "The acquisition of wealth is no longer "
        "the driving force in our lives. "
        "We work to better ourselves and the rest of Humanity."
    ), encoding='utf-8',
)
'''
Arbitrary sequence of bytes to be archived.
'''

# ....................{ TESTS                              }....................
@pytest.mark.parametrize(('filetype',), ARCHIVE_FILETYPES,)
def test_archives_read_write_bytes(
    betse_temp_dir: 'LocalPath', filetype: str) -> None:
    '''
    Unit test both the :func:`reading_bytes` and :func:`writing_bytes` functions of
    the :mod:`betse.util.path.archives` submodule for the passed archive
    filetype guaranteed to be supported by the active Python interpreter.

    Parameters
    ----------
    betse_temp_dir : LocalPath
        Object encapsulating a temporary directory isolated to the current test.
    filetype : str
        Filetype of the archive format to be tested (e.g., "bz2", "gz", "zip").
    '''

    # Defer heavyweight imports.
    from betse.util.path import archives

    # Absolute path of an archive file of arbitrary basename and the passed
    # filetype in this temporary directory.
    archive_filepath = betse_temp_dir.join('memory_alpha.' + filetype)

    # This path as a string rather than "LocalPath" object.
    archive_filename = str(archive_filepath)

    # Ensure this archive filetype to be supported.
    assert archives.is_filetype(archive_filename)

    # Create this archive.
    with archives.write_bytes(archive_filename) as archive_writer:
        archive_writer.write(ARCHIVE_BYTES)

    # Assert this archive to have been created.
    assert archive_filepath.check(file=1)

    # Unpack this archive.
    with archives.read_bytes(archive_filename) as archive_reader:
        archive_bytes_read = archive_reader.read()

    # Assert this archive to contain the bytes written to that archive.
    assert ARCHIVE_BYTES == archive_bytes_read
