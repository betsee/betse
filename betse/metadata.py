#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2025 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Metadata constants synopsizing high-level application behaviour.

Design
------
Metadata constants defined by this submodule are intentionally *not* defined as
metadata properties of the :class:`betse.util.app.meta.appmetaabc` abstract base
class. Why? Because doing so would prevent their use from the top-level
``setup.py`` scripts defined by downstream consumers (e.g., BETSEE GUI), which
would render these constants effectively useless for their principal use case.
'''

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# CAUTION: This submodule has largely (but *NOT* entirely) been obsoleted by the
# top-level "pyproject.toml" file, which should be strongly preferred and
# towards which we should refactor this project away from this submodule.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# ....................{ IMPORTS                            }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To avoid race conditions during setuptools-based installation, this
# module may import *ONLY* from modules guaranteed to exist at the start of
# installation. This includes all standard Python and application modules but
# *NOT* third-party dependencies, which if currently uninstalled will only be
# installed at some later time in the installation.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

import sys

# ....................{ METADATA                           }....................
NAME = 'BETSE'
'''
Human-readable application name.
'''


LICENSE = '2-clause BSD'
'''
Human-readable name of the license this application is licensed under.
'''

# ....................{ PYTHON ~ version                   }....................
PYTHON_VERSION_MIN = '3.11.0'
'''
Human-readable minimum version of Python required by this application as a
``.``-delimited string.
'''


def _convert_version_str_to_tuple(version_str: str) -> tuple:
    '''
    Convert the passed human-readable ``.``-delimited version string into a
    machine-readable version tuple of corresponding integers.
    '''
    assert isinstance(version_str, str), (
        '"{}" not a version string.'.format(version_str))

    return tuple(int(version_part) for version_part in version_str.split('.'))


PYTHON_VERSION_MIN_PARTS = _convert_version_str_to_tuple(PYTHON_VERSION_MIN)
'''
Machine-readable minimum version of Python required by this application as a
tuple of integers.
'''


# Validate the version of the active Python interpreter *BEFORE* subsequent
# code possibly depending on this version. Since this version should be
# validated both at setuptools-based install time and post-install runtime
# *AND* since this module is imported sufficiently early by both, stash this
# validation here to avoid duplication of this logic and hence the hardcoded
# Python version.
#
# The "sys" module exposes three version-related constants for this purpose:
#
# * "hexversion", an integer intended to be specified in an obscure (albeit
#   both efficient and dependable) hexadecimal format: e.g.,
#    >>> sys.hexversion
#    33883376
#    >>> '%x' % sys.hexversion
#    '20504f0'
# * "version", a human-readable string: e.g.,
#    >>> sys.version
#    2.5.2 (r252:60911, Jul 31 2008, 17:28:52)
#    [GCC 4.2.3 (Ubuntu 4.2.3-2ubuntu7)]
# * "version_info", a tuple of three or more integers *OR* strings: e.g.,
#    >>> sys.version_info
#    (2, 5, 2, 'final', 0)
#
# For sanity, this application will *NEVER* conditionally depend upon the
# string-formatted release type of the current Python version exposed via the
# fourth element of the "version_info" tuple. Since the first three elements of
# that tuple are guaranteed to be integers *AND* since a comparable 3-tuple of
# integers is declared above, comparing the former and latter yield the
# simplest and most reliable Python version test.
#
# Note that the nearly decade-old and officially accepted PEP 345 proposed a
# new field "requires_python" configured via a key-value pair passed to the
# call to setup() in "setup.py" (e.g., "requires_python = ['>=2.2.1'],"), that
# field has yet to be integrated into either disutils or setuputils. Hence,
# that field is validated manually in the typical way.
if sys.version_info[:3] < PYTHON_VERSION_MIN_PARTS:
    # Human-readable current version of Python. Ideally, "sys.version" would be
    # leveraged here instead; sadly, that string embeds significantly more than
    # merely a version and hence is inapplicable for real-world usage: e.g.,
    #
    #     >>> import sys
    #     >>> sys.version
    #     '3.6.5 (default, Oct 28 2018, 19:51:39) \n[GCC 7.3.0]'
    PYTHON_VERSION = '.'.join(
        str(version_part) for version_part in sys.version_info[:3])

    # Die ignominiously.
    raise RuntimeError(
        '{} requires at least Python {}, but the active interpreter '
        'is only Python {}. We feel deep sadness for you.'.format(
            NAME, PYTHON_VERSION_MIN, PYTHON_VERSION))

# ....................{ METADATA ~ version                 }....................
VERSION = '1.5.1'
'''
Human-readable application version as a ``.``-delimited string.
'''


CODENAME = 'Nicest Nestor'
'''
Human-readable code name associated with the current version of this
application.

This code name consists of an arbitrary adjective followed by the last name of
an arbitrary academic associated with field of bioelectricity whose first
letter is the same as the first letter of that adjective.

See Also
----------
:download:`/doc/rst/RELEASE.md`
    Further details on this incredible code name schema.
'''


VERSION_PARTS = _convert_version_str_to_tuple(VERSION)
'''
Machine-readable application version as a tuple of integers.
'''

# ....................{ METADATA ~ tag                     }....................
GIT_TAG_COMPAT_OLDEST = 'v0.5.0'
'''
Git-specific tag of the oldest version of this application for which the
current version of this application guarantees backward compatibility.

In this context, the term "backward compatibility" means the capacity for all
versions of this application newer than this tag to successfully load all:

* Simulation configuration files loadable by this older version.
* Pickled seeds, initializations, and simulations saved by this older version.

This specific version was selected for the simple reason that this was the
first version of this application to guarantee backward compatibility, as
implemented by the :mod:`betse.science.compat.compatconf` submodule.

See Also
----------
:func:`betse_test.func.sim.test_sim.test_cli_sim_compat`
    Functional test in our test suite programmatically guaranteeing backward
    compatibility with this older version.
'''

# ....................{ METADATA ~ synopsis                }....................
# Note that a human-readable multiline description is exposed via the top-level
# "setup.py" script. This description is inefficiently culled from the contents
# of the top-level "README.rst" file and hence omitted here. (Doing so here
# would significantly increase program startup costs with little to no benefit.)
SYNOPSIS = 'BETSE, the BioElectric Tissue Simulation Engine.'
'''
Human-readable single-line synopsis of this application.

By PyPI design, this string must *not* span multiple lines or paragraphs.
'''


DESCRIPTION = (
    'The BioElectric Tissue Simulation Engine (BETSE) is a '
    'discrete exterior calculus simulator for '
    '2D computational multiphysics problems in '
    'the life sciences -- including '
    '(electro)diffusion, '
    '(electro)osmosis, '
    'galvanotaxis, '
    'voltage-gated ion channels, '
    'gene regulatory networks, '
    'and biochemical reaction networks.'
)
'''
Human-readable multiline description of this application.

By :mod:`argparse` design, this string may (and typically should) span both
multiple lines and paragraphs. Note that this string is *not* published to
PyPI, which accepts reStructuredText (rst) and is thus passed the contents of
the top-level :doc:`/README` file instead.
'''

# ....................{ METADATA ~ authors                 }....................
AUTHORS = 'Alexis Pietak, Cecil Curry, et al.'
'''
Human-readable list of all principal authors of this application as a
comma-delimited string.

For brevity, this string *only* lists authors explicitly assigned copyrights.
For the list of all contributors regardless of copyright assignment or
attribution, see the top-level `AUTHORS.md` file.
'''


AUTHOR_EMAIL = 'alexis.pietak@gmail.com'
'''
Email address of the principal corresponding author (i.e., the principal author
responding to public correspondence).
'''

# ....................{ METADATA ~ urls                    }....................
URL_HOMEPAGE = 'https://github.com/betsee/betse'
'''
URL of this application's homepage.
'''


URL_DOWNLOAD = f'{URL_HOMEPAGE}/archive/{VERSION}.tar.gz'
'''
URL of the source tarball for the current version of this application.

This URL assumes a tag whose name is ``v{VERSION}`` where ``{VERSION}`` is the
human-readable current version of this application (e.g., ``v0.4.0``) to exist.
Typically, no such tag exists for live versions of this application -- which
have yet to be stabilized and hence tagged. Hence, this URL is typically valid
*only* for previously released (rather than live) versions of this application.
'''


URL_ISSUES = f'{URL_HOMEPAGE}/issues'
'''
URL of this package's issue tracker.
'''


URL_RELEASES = f'{URL_HOMEPAGE}/releases'
'''
URL of this package's release list.
'''

# ....................{ METADATA ~ python                  }....................
PACKAGE_NAME = NAME.lower()
'''
Fully-qualified name of the top-level Python package implementing this
application.

Caveats
----------
**Prefer the application-agnostic
:meth:`betse.appmeta.app_meta.package_name` property instead.** This
application-specific global fails to generalize to downstream consumers (e.g.,
BETSEE) and hence is usable *ONLY* for low-level installation-time use cases.
'''
