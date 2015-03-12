#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Metadata constants synopsizing high-level `betse` behaviour.

Unrelatedly, this module also validates the version of the active Python 3
interpreter. An exception is raised if such version is insufficient.
'''

# ....................{ IMPORTS                            }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To avoid race conditions during setuptools-based installation, this
# module may import *ONLY* from packages guaranteed to exist at installation
# time (e.g., stock Python packages).
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# ....................{ METADATA                           }....................
# General-purpose metadata.

NAME = 'BETSE'
'''Human-readable program name.'''

DESCRIPTION = ''.join((
    NAME, ', the [B]io[E]lectric [T]issue [S]imulation [E]ngine: '
    'bio-realistic modelling of dynamic electrochemical phenomena in '
    'gap junction-networked cell collectives, ',
    'with a focus on spatio-temporal pattern formation.',
))
'''Human-readable program description.'''

AUTHORS = 'Alexis Pietak, Cecil Curry'
# AUTHORS = 'Alexis Pietak, Cecil Curry, et al.'
'''Human-readable program authors as a comma-delimited list.'''

# ....................{ METADATA ~ scripts                 }....................
SCRIPT_NAME_CLI = NAME.lower()
'''
Basename of the CLI-specific Python script wrapper created by `setuptools`
installation.
'''

SCRIPT_NAME_GUI = SCRIPT_NAME_CLI + '-qt'
'''
Basename of the GUI-specific Python script wrapper created by `setuptools`
installation.
'''

# ....................{ METADATA ~ versions                }....................
__version__ = '0.1'
'''
Version specifier.

For PEP 8 compliance, such specifier is exposed as the canonical variable
`__variable__` rather than a typical constant (e.g., `VARIABLE`).
'''

# Program version as a tuple adhering to this semi-standard version specifier.
__version_info__ = tuple(
    int(version_part) for version_part in __version__.split('.'))
'''
Version specifier as a tuple of integers rather than `.`-delimited string.

For PEP 8 compliance, such specifier is exposed as the canonical variable
`__variable_info__` rather than a typical constant (e.g., `VARIABLE_PARTS`).
'''

# ....................{ METADATA ~ setuptools              }....................
# setuptools-specific metadata required outside of setuptools-based
# installations, typically for performing runtime validation of the current
# Python environment.

DEPENDENCY_SETUPTOOLS = 'setuptools >= 5.0'
'''
Version of `setuptools` required by `betse` at both install and runtime.

For simplicity, such version is a `setuptools`-specific requirements string.
'''

#FIXME: Reenable commented dependencies when actually imported by the codebase.

# Such dependencies are also checked for importability as modules at runtime and
# hence case-sensitive (e.g., "PySide" is importable but "pyside" is *NOT*).
DEPENDENCIES_RUNTIME = [
    # setuptools is currently required at both install and runtime. At runtime,
    # setuptools is used to validate that required dependencies are available.
    DEPENDENCY_SETUPTOOLS,

    # matplotlib 1.4.0 switched Python compatability layers from "2to3" to
    # "six". PyInstaller does *NOT* yet support the latter, preventing freezing
    # and hence end-user use of such version of matplotlib.
    'matplotlib >= 1.3.0',
    # 'matplotlib < 1.4.0',

    # Dependencies directly required by BETSE.
    # 'PySide >= 1.2.0',
    'numpy >= 1.8.0',
    'scipy >= 0.12.0',
    # 'voluptuous >= 0.8.7',
    'yaml >= 3.10',

    # Dependencies transitively but *NOT* directly required by BETSE. To detect
    # missing such dependencies in a human-readable manner, such dependencies
    # are explicitly listed as well.
    'six >= 1.8.0',
]
'''
Set of all mandatory runtime dependencies for `betse`.

For simplicity, such set is formatted as a list of `setuptools`-specific
requirements strings.

See Also
----------
README.md
    Human-readable list of such dependencies.
'''

# ....................{ PYTHON ~ version                   }....................
# Validate the version of the active Python interpreter *BEFORE* subsequent code
# possibly depending on such version. Since such version should be validated
# both at setuptools-based install time and post-install runtime *AND* since
# this module is imported sufficiently early by both, stash such validation here
# to avoid duplication of such logic and hence the hardcoded Python version.
#
# The "sys" module exposes three version-related constants for such purpose:
#
# * "hexversion", an integer intended to be used in hexadecimal format: e.g.,
#    >>> sys.hexversion
#    33883376
#    >>> '%x' % sys.hexversion
#    '20504f0'
# * "version", a human-readable string: e.g.,
#    >>> sys.version
#    2.5.2 (r252:60911, Jul 31 2008, 17:28:52)
#    [GCC 4.2.3 (Ubuntu 4.2.3-2ubuntu7)]
# * "version_info", a tuple of three or more integers or strings: e.g.,
#    >>> sys.version_info
#    (2, 5, 2, 'final', 0)
#
# Only the first lends itself to reliable comparison.
#
# Note that the nearly decade-old and officially accepted PEP 345 proposed a new
# field "requires_python" configured via a key-value pair passed to the call to
# setup() in "setup.py" (e.g., "requires_python = ['>=2.2.1'],"), such field has
# yet to be integrated into either disutils or setuputils. Hence, such field is
# validated manually in the typical way. Behead the infidel setuptools!
import sys
if sys.hexversion < 0x03030000:
    raise RuntimeError(''.join((
        NAME, ' ',
        'requires at least Python 3.3. However, the active Python ',
        'interpreter is of version:\n\n',
        sys.version, '\n\n',
        'We feel sadness for you.',
    )))

# --------------------( WASTELANDS                         )--------------------
#FUXME: There appears to be no way to provide build-time dependencies to
#setuptools. Nonetheless, there should be, so we preserve this for posterity.

# DEPENDENCIES_INSTALL = [
#     'setuptools >= 7.0',
# ]
# '''
# Set of all mandatory install-time dependencies for `betse`.
#
# For simplicity, such set is formatted as a list of `setuptools`-specific
# requirements strings.
#
# See Also
# ----------
# README.md
#     Human-readable list of such dependencies.
# '''

# See "README.md" for human-readable versions of such dependency lists.
# For simplicity, such set is formatted as a list of `setuptools`-specific
# requirements strings.

    #FUXME: Reenable after installation.
# COMMAND_NAME_PREFIX = NAME.lower()
# '''
# Substring prefixing the basenames of all `betse`-specific Python script wrappers
# created by `setuptools` installation.
# '''

#  In particular, this implies that
# packages possibly installed by such installation (e.g., mandatory
# dependencies) must *NOT* be imported from.
#
# Yes, this is terrible. No, we didn't make the ad-hoc rules.

# Such type is validated by comparison of
# And yes: this *IS* the typical way
# for validating Python versions. Note that
