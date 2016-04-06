#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Metadata constants synopsizing high-level `betse` behaviour.

Python Version
----------
For uniformity between both the main BETSE codebase and the "setup.py"
setuptools script importing this module, this module also validates the version
of the active Python 3 interpreter. An exception is raised if such version is
insufficient.

BETSE currently requires **Python 3.4**, as:

* Python 3.3 provides insufficient machinery for dynamically inspecting modules
  at runtime. In particular, both the long-standing `imp.find_module()` function
  and the `importlib.find_loader()` function introduced by Python 3.3 require
  all parent packages of the passed module to be recursively imported _before_
  such functions are called; failing to do so results in such functions
  unconditionally returning `None`. Since this has been the source of numerous
  subtle issues throughout the codebase, Python 3.3 is strictly out. Since
  modern Linux distributions have long since switched away from Python 3.3 as
  their default Python 3 interpreters, this _should_ impose no hardships.
'''

# ....................{ IMPORTS                            }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To avoid race conditions during setuptools-based installation, this
# module may import *ONLY* from packages guaranteed to exist at the start of
# installation. This includes all standard Python and BETSE packages but *NOT*
# third-party dependencies, which if currently uninstalled will only be
# installed at some later time in the installation.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

import sys

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
#
# WARNING: When modifying this, please document the justification for doing so
# in the "Python Version" subsection of this module's docstring above.
if sys.hexversion < 0x03040000:
    raise RuntimeError(''.join((
        NAME, ' ',
        'requires at least Python 3.4, ',
        'but the active Python interpreter is only\n',
        'Python ', sys.version, '. We feel sadness for you.',
    )))

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
__version__ = '0.4.0'
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

DEPENDENCY_SETUPTOOLS = 'setuptools >= 3.3'
'''
Version of `setuptools` required by `betse` at both install and runtime.

For simplicity, such version is a `setuptools`-specific requirements string.

Note that `betse` only requires the `pkg_resources` package installed along
with `setuptools` rather than the `setuptools` package itself. Since there
appears to exist no means of asserting a dependency on only `pkg_resources`,
we pretend to require `setuptools` itself. This is non-ideal, of course.
'''

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# NOTE: Changes to dependency names declared by this set *MUST* be synchronized
# with the corresponding keys of the
# "betse.util.py.modules.SETUPTOOLS_PROJECT_TO_MODULE_NAME" dictionary.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
DEPENDENCIES_RUNTIME = [
    # setuptools is currently required at both install and runtime. At runtime,
    # setuptools is used to validate that required dependencies are available.
    DEPENDENCY_SETUPTOOLS,

    # Dependencies directly required by BETSE.
    'Matplotlib >= 1.4.0',
    'Numpy >= 1.8.0',
    'SciPy >= 0.12.0',
    'PyYAML >= 3.10',
    # 'PySide >= 1.2.0',
    # 'Yamale >= 1.5.0',

    # Dependencies transitively but *NOT* directly required by BETSE. To detect
    # missing dependencies in a human-readable manner, these dependencies are
    # explicitly listed as well.
    'six >= 1.5.2',      # required by everything that should not be
    'Pillow >= 2.3.0',   # required by the "scipy.misc.imread" module
]
'''
Set of all mandatory runtime dependencies for `betse`.

For simplicity, this set is formatted as a list of `setuptools`-specific
requirements strings whose:

* First word is the name of the `setuptools`-specific project being required,
  which may have no relation to the name of that project's top-level module or
  package (e.g., the `PyYAML` project's top-level package is `yaml`). For human
  readability in error messages, this name should typically be case-sensitively
  capitalized -- despite being parsed case-insensitively by `setuptools`.
* Second word is a numeric comparison operator.
* Third word is the version specifier of that project required by that
  comparison.

See Also
----------
README.md
    Human-readable list of these dependencies.
'''
