#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Metadata constants synopsizing high-level application behaviour.

Python Version
----------
For uniformity between this codebase and the ``setup.py`` setuptools script
importing this module, this module also validates the version of the active
Python 3 interpreter. An exception is raised if this version is insufficient.

This application currently requires **Python 3.5.** This application previously
required only Python 3.4.0 for access to the standard :mod:`enum` module
introduced by that major release. Since Python 3.4.0 is approaching its
official end-of-life (EOL) and scheduled to be decomissioned shortly, this
application now unavoidably requires the next major release of Python despite
currently leveraging no features introduced by that release.

Note that Python < 3.4.0:

* Provides insufficient machinery for dynamically inspecting modules at
  runtime. In particular, both the long-standing :func:`imp.find_module`
  function and the :func:`importlib.find_loader` function introduced by Python
  3.3 require all parent packages of the passed module to be recursively
  imported *before* these functions are called; failing to do so results in
  these functions unconditionally returning ``None``. Since this has been the
  source of numerous subtle issues throughout this codebase, Python 3.3 is
  strictly out. Since most modern Linux distributions have adopted Python 3.4
  as the default Python
  3 interpreters, this *should* impose no hardship.
* Fails to provide the :mod:`enum` module introduced by Python 3.4, which both
  standardizes and simplifies enumeration implementations.

Design
----------
Metadata constants defined by this submodule are intentionally *not* defined as
metadata properties of the :class:`betse.util.app.meta.appmetaabc` abstract base
class. Why? Because doing so would prevent their use from the top-level
``setup.py`` scripts defined by downstream consumers (e.g., BETSEE GUI), which
would render these constants effectively useless for their principal use case.
'''

# ....................{ IMPORTS                           }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To avoid race conditions during setuptools-based installation, this
# module may import *ONLY* from modules guaranteed to exist at the start of
# installation. This includes all standard Python and application modules but
# *NOT* third-party dependencies, which if currently uninstalled will only be
# installed at some later time in the installation.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

import sys

# ....................{ METADATA                          }....................
NAME = 'BETSE'
'''
Human-readable application name.
'''


LICENSE = '2-clause BSD'
'''
Human-readable name of the license this application is licensed under.
'''

# ....................{ PYTHON ~ version                  }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: Changes to this section *MUST* be synchronized with:
# * The corresponding section of the "betsee.metadata" submodule.
# * Front-facing documentation (e.g., "README.rst", "doc/md/INSTALL.md").
# On bumping the minimum required version of Python, consider also documenting
# the justification for doing so in the "Python Version" section of this
# submodule's docstring above.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PYTHON_VERSION_MIN = '3.5.0'
'''
Human-readable minimum version of Python required by this application as a
``.``-delimited string.

See Also
----------
"Python Version" section of this submodule's docstring for a detailed
justification of this constant's current value.
'''


PYTHON_VERSION_MINOR_MAX = 8
'''
Maximum minor stable version of this major version of Python currently released
(e.g., ``5`` if Python 3.5 is the most recent stable version of Python 3.x).
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

# ....................{ METADATA ~ version                }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: When modifying the current version of this application below,
# consider adhering to the Semantic Versioning schema. Specifically, the
# version should consist of three "."-delimited integers
# "{major}.{minor}.{patch}", where:
#
# * "{major}" specifies the major version, incremented only when either:
#   * Breaking configuration file backward compatibility. Since this
#     application's public API is its configuration file format rather than a
#     subset of the code itself (e.g., public submodules, classes), no change
#     to the code itself can be considered to break backward compatibility
#     unless that change breaks the configuration file format.
#   * Implementing headline-worthy functionality (e.g., a GUI). Technically,
#     this condition breaks the Semantic Versioning schema, which stipulates
#     that *ONLY* changes breaking backward compatibility warrant major bumps.
#     But this is the real world. In the real world, significant improvements
#     are rewarded with significant version changes.
#   In either case, the minor and patch versions both reset to 0.
# * "{minor}" specifies the minor version, incremented only when implementing
#   customary functionality in a manner preserving such compatibility. In this
#   case, the patch version resets to 0.
# * "{patch}" specifies the patch version, incremented only when correcting
#   outstanding issues in a manner preserving such compatibility.
#
# When in doubt, increment only the minor version and reset the patch version.
# For further details, see http://semver.org.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

VERSION = '1.1.1'
'''
Human-readable application version as a ``.``-delimited string.
'''


CODENAME = 'Nicer Nestor'
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

# ....................{ METADATA ~ tag                    }....................
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

# ....................{ METADATA ~ synopsis               }....................
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

# ....................{ METADATA ~ authors                }....................
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

# ....................{ METADATA ~ urls                   }....................
URL_HOMEPAGE = 'https://gitlab.com/betse/betse'
'''
URL of this application's homepage.
'''


URL_DOWNLOAD = '{}/repository/archive.tar.gz?ref=v{}'.format(
    URL_HOMEPAGE, VERSION)
'''
URL of the source tarball for the current version of this application.

This URL assumes a tag whose name is ``v{VERSION}`` where ``{VERSION}`` is the
human-readable current version of this application (e.g., ``v0.4.0``) to exist.
Typically, no such tag exists for live versions of this application -- which
have yet to be stabilized and hence tagged. Hence, this URL is typically valid
*only* for previously released (rather than live) versions of this application.
'''

# ....................{ METADATA ~ python                 }....................
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
