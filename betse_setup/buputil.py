#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Utility and convenience functions for application-specific :mod:`setuptools`
subcommands.

Design
----------
This method intentionally duplicates existing utility functions provided by the
:mod:`betse.util` subpackage. While duplication breaks DRY ("Don't Repeat
Yourself") and hence is usually harmful, there are valid reasons to do so here.
Namely, :mod:`betse.util` functionality:

* Assumes BETSE to be available. While this is certainly the case when this
  file resides in the BETSE codebase, this is *not* necessarily the case when
  this file is copied into and hence resides in the codebases of other projects
  (e.g., BETSEE). In these projects, BETSE is merely yet another dependency
  that is typically unavailable at installation time.
* Raises BETSE-specific exceptions rooted at the BETSE-specific
  :class:`betse.exception.BetseException` superclass. :mod:`setuptools`
  subcommands, on the other hand, are expected to only raise
  :mod:`distutils`-specific exceptions rooted at the :mod:`distutils`-specific
  :class:`DistutilsError` superclass.
* Assumes logging to be configured. :mod:`setuptools`, however, assumes
  logging to *not* be configured -- and provides no assistance in doing so.
* Could theoretically import third-party dependencies unavailable at
  :mod:`setuptools subcommand time (e.g., due to the ``install`` or ``develop``
  subcommands *not* having been run yet). While no :mod:`betse.util` submodules
  should do so, the horrid possibility remains.

Since duplicating these functions here is no significant maintenance burden
*and* since attempting to reuse these functions here would introduce spurious
side effects, we adopt the former approach.
'''

# ....................{ IMPORTS                           }....................
import setuptools, sys
from distutils.version import StrictVersion

# ....................{ EXCEPTIONS                        }....................
def die_unless_setuptools_version_at_least(
    setuptools_version_min: str) -> None:
    '''
    Raise an exception unless the currently installed version of
    :mod:`setuptools` is at least as recent as the passed minimum version.

    Parameters
    ----------
    setuptools_version_min : str
        Human-readable ``.``-delimited specifier of the minimum version of
        :mod:`setuptools` required at installation time by this application.

    Raises
    ----------
    Exception
        If the currently installed version of :mod:`setuptools` is older than
        the passed minimum version.
    '''
    assert isinstance(setuptools_version_min, str), (
        '"{}" not a string.'.format(setuptools_version_min))

    # If the currently installed version of setuptools is older than this
    # minimum version, raise an exception.
    if (
        StrictVersion(setuptools.__version__) <
        StrictVersion(setuptools_version_min)
    ):
        raise Exception(
            'setuptools >= {} required by this application, but only '
            'setuptools {} found.'.format(
                setuptools_version_min, setuptools.__version__))

# ....................{ GETTERS                           }....................
def get_chars(filename: str, encoding: str = 'utf-8') -> str:
    '''
    String of all characters contained in the plaintext file with the passed
    filename encoded with the passed encoding.

    Parameters
    ----------
    filename : str
        Relative or absolute path of the plaintext text to be read.
    encoding : optional[str]
        Name of the encoding to be used. Defaults to UTF-8.

    Returns
    ----------
    str
        String of all characters decoded from this file's byte content.
    '''
    assert isinstance(filename, str), '"{}" not a string.'.format(filename)
    assert isinstance(encoding, str), '"{}" not a string.'.format(encoding)

    with open(filename, mode='rt', encoding=encoding) as text_file:
        return text_file.read()


def get_description() -> str:
    '''
    Human-readable multiline description of this application in
    reStructuredText (reST) format.

    To minimize synchronization woes, this description is identical to the
    contents of the :doc:`/README.rst` file. When submitting this application
    package to PyPI, this description is re-used verbatim as this package's
    front matter.

    Caveats
    ----------
    This function is I/O intensive and hence should be called sparingly --
    ideally, only once by this application's top-level ``setup.py`` script.
    '''

    # Relative path of this application's front-facing documentation in
    # reStructuredText format, required by PyPI. This path resides outside this
    # application's package tree and hence is inlined here rather than provided
    # by the "betsee.guiappmeta" submodule.
    DESCRIPTION_FILENAME = 'README.rst'

    # Description read from this description file.
    try:
        description = get_chars(DESCRIPTION_FILENAME)
        # print('description: {}'.format(_DESCRIPTION))
    # If this file is *NOT* readable, print a non-fatal warning and reduce this
    # description to the empty string. While unfortunate, this description is
    # *NOT* required for most operations and hence mostly ignorable.
    except Exception as exception:
        description = ''
        _output_warning(
            'Description file "{}" not found or unreadable:\n{}'.format(
                DESCRIPTION_FILENAME, exception))

    # Retcurn this description.
    return description

# ....................{ SANITIZERS                        }....................
def sanitize_classifiers(
    classifiers: list,
    python_version_min_parts: tuple,
    python_version_minor_max: int,
) -> list:
    '''
    List of all PyPI-specific trove classifier strings synopsizing this
    application, manufactured by appending classifiers synopsizing this
    application's support for Python major versions (e.g.,
    ``Programming Language :: Python :: 3.6``, a classifier implying this
    application to successfully run under Python 3.6) to the passed list.

    Parameters
    ----------
    classifiers : list
        List of all PyPI-specific trove classifier strings to be sanitized.
    python_version_min_parts : tuple
        Minimum fully-specified version of Python required by this application
        as a tuple of integers (e.g., ``(3, 5, 0)`` if this application
        requires at least Python 3.5.0).
    python_version_minor_max : int
        Maximum minor stable version of the current Python 3.x mainline (e.g.,
        ``9`` if Python 3.9 is the most recent stable version of Python 3.x).

    Returns
    ----------
    list
        List of all sanitized PyPI-specific trove classifier strings.
    '''
    assert isinstance(classifiers, list), '"{}" not a list.'.format(
        classifiers)
    assert isinstance(python_version_min_parts, tuple), (
        '"{}" not a tuple.'.format(python_version_min_parts))
    assert isinstance(python_version_minor_max, int), (
        '"{}" not an integer.'.format(python_version_minor_max))

    # Major version of Python required by this application.
    PYTHON_VERSION_MAJOR = python_version_min_parts[0]

    # List of classifiers to return, copied from the passed list for safety.
    classifiers_sane = classifiers[:]

    # For each minor version of Python 3.x supported by this application,
    # formally classify this version as such.
    for python_version_minor in range(
        python_version_min_parts[1], python_version_minor_max):
        classifiers.append(
            'Programming Language :: Python :: {}.{}'.format(
                PYTHON_VERSION_MAJOR, python_version_minor,))
    # print('classifiers: {}'.format(_CLASSIFIERS))

    # Return this sanitized list of classifiers.
    return classifiers_sane

# ....................{ OUTPUTTERS                        }....................
def _output_warning(*warnings) -> None:
    '''
    Print the passed warning messages to standard error.
    '''

    print('WARNING: ', *warnings, file=sys.stderr)
