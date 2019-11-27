#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Metadata constants synopsizing high-level application dependencies.

Design
----------
Metadata constants defined by this submodule are intentionally *not* defined as
metadata properties of the :class:`betse.util.app.meta.appmetaabc` abstract
base class. Why? Because doing so would prevent their use from the top-level
``setup.py`` scripts defined by downstream consumers (e.g., BETSEE), which
would render these constants effectively useless for their principal use case.
'''

#FIXME: Setuptools should no longer be required at runtime under Python >= 3.8.
#Why? Because Python 3.8 wisely introduced the new "importlib.metadata" module,
#enabling inspection of package manager-installed metadata (e.g., version,
#dependencies, entry points) *WITHOUT* requiring that manager to be installed.
#"importlib.metadata" is likely to be substantially more efficient than the
#setuptools-specific "pkg_resources" and "setuptools" packages, which are
#well-known to be CPU bottlenecks. Of course, supporting this sanely will
#necessitate that we refactor the codebase as follows:
#
#* Define a new "betse.util.py.pyproject" submodule, providing high-level
#  access to project-specific metadata rather than lower-level module- or
#  package-specific metadata. This submodule should define functions generally
#  resembling those already defined by the "betse.lib.setuptools.setuptool"
#  submodule as follows:
#  * If this is Python >= 3.8, these functions should be implemented in terms
#    of the new "importlib.metadata" submodule.
#  * If this is Python < 3.8, these functions should be implemented in terms
#    of the old "betse.lib.setuptools.setuptool" submodule.
#
#For backward compatibility, the "betse.lib.setuptools.setuptool" submodule
#should be preserved in perpetuity as is. It works, so don't break it! Please.

# ....................{ IMPORTS                           }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To avoid race conditions during setuptools-based installation, this
# module may import *ONLY* from modules guaranteed to exist at the start of
# installation. This includes all standard Python and application modules but
# *NOT* third-party dependencies, which if currently uninstalled will only be
# installed at some later time in the installation.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

from collections import namedtuple

# ....................{ LIBS ~ install : mandatory        }....................
#FIXME: The current situation of "setuptools" versioning duplication detailed
#below is absolutely insane and incentive enough for us to gradually transition
#away from the "setuptools"-based "setup.py" ecosystem to the "poetry"-based
#"pyproject.toml" ecosystem. Note, however, that the latter is still largely in
#flux -- which suggests a gestational waiting period. Early adopters always get
#burned hard... and we've been burned hard enough. Even after "poetry"
#solidifies, however, the transition is likely to be non-trivial. Requisite
#changes probably include:
#
#* Embedding the "betse" package under a new top-level "src/" directory.
#* Shifting the contents of both the "betse.metadata.py" *AND* "betse.metadeps"
#  submodules into a "pyproject.toml" file. Note that doing so may
#  fundamentally break BETSEE's current merger algorithm of these files -- an
#  issue that may have no trivial resolution.
#* Refactor both the "betse.metadata.py" *AND* "betse.metadeps" submodules to
#  dynamically deserialize the "pyproject.toml" file and cache the contents of
#  that file into Python-accessible global module variables.

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To avoid dependency conflicts between "pip", "setuptools", BETSE,
# and BETSEE, the value of this global variable *MUST* be synchronized (i.e.,
# copied) across numerous files in both codebases. Specifically, the following
# strings *MUST* be identical:
# * "betse.metadeps.SETUPTOOLS_VERSION_MIN".
# * "betsee.guimetadeps.SETUPTOOLS_VERSION_MIN".
# * The "requires" setting in:
#   * "betse/tox.ini".
#   * "betsee/tox.ini".
# * The "build-backend" setting in:
#   * "betse/pyproject.toml".
#   * "betsee/pyproject.toml".
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# This public global is externally referenced by "setup.py".
SETUPTOOLS_VERSION_MIN = '38.2.0'
'''
Minimum version of :mod:`setuptools` required at both application install- and
runtime as a human-readable ``.``-delimited string.

Motivation
----------
This version derives from both:

* BETSE requirements. BETSE requires the
  :meth:`setuptools.command.easy_install.ScriptWriter.get_args` class method
  and hence at least the oldest version of :mod:`setuptools` to have this
  method. Since official setuptools documentation fails to specify the exact
  version that first defined this method, we fallback to a sufficiently old
  version from 2017 known to define this method: :mod:`setuptools` >= 36.7.2.
* BETSEE requirements. BETSEE requires :mod:`PySide2`, which is distributed as
  a wheel and thus requires wheel support, which in turns requires either
  ``pip`` >= 1.4.0 or :mod:`setuptools` >= 38.2.0. While ``pip`` 1.4.0 is
  relatively ancient, :mod:`setuptools` 38.2.0 is comparatively newer. If the
  current version of :mod:`setuptools` is *not* explicitly validated at
  installation time, older :mod:`setuptools` versions fail on attempting to
  install :mod:`PySide2` with non-human-readable fatal errors resembling:

    $ sudo python3 setup.py develop
    running develop
    running egg_info
    writing betsee.egg-info/PKG-INFO
    writing dependency_links to betsee.egg-info/dependency_links.txt
    writing entry points to betsee.egg-info/entry_points.txt
    writing requirements to betsee.egg-info/requires.txt
    writing top-level names to betsee.egg-info/top_level.txt
    reading manifest template 'MANIFEST.in'
    writing manifest file 'betsee.egg-info/SOURCES.txt'
    running build_ext
    Creating /usr/lib64/python3.6/site-packages/betsee.egg-link (link to .)
    Saving /usr/lib64/python3.6/site-packages/easy-install.pth
    Installing betsee script to /usr/bin
    changing mode of /usr/bin/betsee to 755

    Installed /home/leycec/py/betsee
    Processing dependencies for betsee==0.9.2.0
    Searching for PySide2
    Reading https://pypi.python.org/simple/PySide2/
    No local packages or working download links found for PySide2
    error: Could not find suitable distribution for Requirement.parse('PySide2')

Since these two versions superficially conflict, this string global reduces to
the version satisfying both constraints (i.e., the newest of these two
versions). Equivalently:

.. code-block:: console

   SETUPTOOLS_VERSION_MIN == max('36.7.2', '38.2.0') == '38.2.0'
'''

# ....................{ LIBS ~ runtime : mandatory        }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: Changes to this subsection *MUST* be synchronized with:
# * Front-facing documentation (e.g., "doc/md/INSTALL.md").
# * The "betse.util.py.module.pymodname.DISTUTILS_PROJECT_NAME_TO_MODULE_NAME"
#   dictionary, converting between the setuptools-specific names listed below
#   and the Python-specific module names imported by this application.
# * Gitlab-CI configuration (e.g., the top-level "requirements-conda.txt" file).
# * Third-party platform-specific packages (e.g., Gentoo Linux ebuilds).
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RUNTIME_MANDATORY = {
    # setuptools is currently required at both install and runtime. At runtime,
    # setuptools is used to validate that dependencies are available.
    'setuptools': '>= ' + SETUPTOOLS_VERSION_MIN,

    # Dependencies directly required by this application. Notably:
    #
    # * Numpy 1.13.0 first introduced the optional "axis" keyword argument to
    #   the numpy.unique() function, which this codebase commonly passes.
    'Numpy':  '>= 1.13.0',
    'Pillow': '>= 2.3.0',
    'SciPy':  '>= 0.12.0',
    'dill':   '>= 0.2.3',

    # Matplotlib >= 1.5.0 is required for the newly added "viridis" colormap.
    'matplotlib': '>= 1.5.0',

    # A relatively modern version of "ruamel.yaml" variants is required.
    # Specifically, this application requires:
    #
    # * At least version 0.15.24 or newer of "ruamel.yaml", which resolves a
    #   long-standing parser issue preventing overly complex YAML markup (such
    #   as ours) from being safely roundtripped:
    #   "0.15.24 (2017-08-09):
    #    * (finally) fixed longstanding issue 23 (reported by Antony Sottile),
    #      now handling comment between block mapping key and value correctly"
    # * The new "ruamel.yaml" API first introduced in 0.15.0. While older
    #   versions strictly adhere to the functional PyYAML-compatible API, newer
    #   versions break backward compatibility by entirely supplanting that API
    #   with a modern object-oriented approach. Supporting both isn't worth the
    #   substantial increase in maintenance debt.
    'ruamel.yaml': '>= 0.15.24',

    # Dependencies indirectly required by this application but only optionally
    # required by dependencies directly required by this application. Since the
    # "setup.py" scripts for the latter do *NOT* list these dependencies as
    # mandatory, these dependencies *MUST* be explicitly listed here.

    # Dependencies directly required by dependencies directly required by this
    # application. While these dependencies need *NOT* be explicitly listed
    # here, doing so improves detection of missing dependencies in a
    # human-readable manner.
    'six': '>= 1.5.2',       # required by everything that should not be
}
'''
Dictionary mapping from the :mod:`setuptools`-specific project name of each
mandatory runtime dependency for this application to the suffix of a
:mod:`setuptools`-specific requirements string constraining this dependency.

To simplify subsequent lookup, these dependencies are contained by a dictionary
rather than a simple set or sequence such that each:

* Key is the name of a :mod:`setuptools`-specific project identifying this
  dependency, which may have no relation to the name of that project's
  top-level module or package (e.g., the ``PyYAML`` project's top-level package
  is :mod:`yaml`). For human readability in error messages, this name should
  typically be case-sensitively capitalized -- despite being parsed
  case-insensitively by :mod:`setuptools`.
* Value is either:

  * ``None`` or the empty string, in which case this dependency is
    unconstrained (i.e., any version of this dependency is sufficient).
  * A string of the form ``{comparator} {version}``, where:

    * ``{comparator}`` is a comparison operator (e.g., ``>=``, ``!=``).
    * ``{version}`` is the required version of this project to compare.

Concatenating each such key and value yields a :mod:`setuptools`-specific
requirements string of the form either ``{project_name}`` or ``{project_name}
{comparator} {version}``.

Official :mod:`setuptools` documentation suggests the ``install_requires`` and
``setup_requires`` keys of the ``setup()`` packaging function to accept only
sequences rather than dictionaries of strings. While undocumented, these keys
*do* actually accept both sequences and dictionaries of strings.

Caveats
----------
This application requires :mod:`setuptools` at both installation time *and*
runtime -- in the latter case, to validate all application dependencies at
runtime. Note that doing so technically only requires the :mod:`pkg_resources`
package installed with :mod:`setuptools` rather than the :mod:`setuptools`
package itself. Since there exists no means of asserting a dependency on only
:mod:`pkg_resources`, however, :mod:`setuptools` is depended upon instead.

See Also
----------
https://setuptools.readthedocs.io/en/latest/setuptools.html#id12
    Further details on the :mod:`setuptools` string format for dependencies.
:download:`/doc/md/INSTALL.md`
    Human-readable list of these dependencies.
'''

# ....................{ LIBS ~ runtime : optional         }....................
RUNTIME_OPTIONAL = {
    # To simplify subsequent lookup at runtime, project names for optional
    # dependencies should be *STRICTLY LOWERCASE*. Since setuptools parses
    # project names case-insensitively, case is only of internal relevance.

    # Optional dependencies directly required by this application.
    'distro':   '>= 1.0.4',
    'pympler':  '>= 0.4.1',
    'ptpython': '>= 0.29',

    #FIXME: Uncomment once eventually used, which is probably inevitable now.
    # 'psutil':   '>= 5.3.0',

    # A relatively modern version of NetworkX *EXCLUDING* 1.11, which
    # critically broke backwards compatibility by coercing use of the
    # unofficial inactive "pydotplus" PyDot fork rather than the official
    # active "pydot" PyDot project, is directly required by this application.
    # NetworkX >= 1.12 reverted to supporting "pydot", thus warranting
    # blacklisting of only NetworkX 1.11. (It is confusing, maybe? Yes!)
    'networkx': '>= 1.8, != 1.11',
    'pydot': '>= 1.0.28',
}
'''
Dictionary mapping from the :mod:`setuptools`-specific project name of each
optional runtime dependency for this application to the suffix of a
:mod:`setuptools`-specific requirements string constraining this dependency.

See Also
----------
:data:`RUNTIME_MANDATORY`
    Further details on dictionary structure.
:download:`/doc/md/INSTALL.md`
    Human-readable list of these dependencies.
:func:`get_dependencies_runtime_optional_tuple`
    Function converting this dictionary of key-value string pairs into a tuple
    of strings (e.g., within :download:`/setup.py`).
'''

# ....................{ LIBS ~ testing : mandatory        }....................
TESTING_MANDATORY = {
    # setuptools is currently required at testing time as well. If ommitted,
    # "tox" commonly fails at venv creation time with exceptions resembling:
    #
    #     GLOB sdist-make: /home/leycec/py/betse/setup.py
    #     py36 inst-nodeps: /home/leycec/py/betse/.tox/.tmp/package/1/betse-1.1.1.zip
    #     ERROR: invocation failed (exit code 1), logfile: /home/leycec/py/betse/.tox/py36/log/py36-3.log
    #     =================================================== log start ===================================================
    #     Processing ./.tox/.tmp/package/1/betse-1.1.1.zip
    #         Complete output from command python setup.py egg_info:
    #         Traceback (most recent call last):
    #           File "<string>", line 1, in <module>
    #           File "/tmp/pip-0j3y5x58-build/setup.py", line 158, in <module>
    #             buputil.die_unless_setuptools_version_at_least(metadeps.SETUPTOOLS_VERSION_MIN)
    #           File "/tmp/pip-0j3y5x58-build/betse_setup/buputil.py", line 74, in die_unless_setuptools_version_at_least
    #             setuptools_version_min, setuptools.__version__))
    #         Exception: setuptools >= 38.2.0 required by this application, but only setuptools 28.8.0 found.
    #
    #         ----------------------------------------
    #     Command "python setup.py egg_info" failed with error code 1 in /tmp/pip-0j3y5x58-build/
    #     You are using pip version 9.0.1, however version 19.3.1 is available.
    #     You should consider upgrading via the 'pip install --upgrade pip' command.
    #
    #     ==================================================== log end ====================================================
    'setuptools': '>= ' + SETUPTOOLS_VERSION_MIN,

    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # WARNING: This py.test requirement *MUST* be manually synchronized to the
    # same requirement in the downstream "betsee.guimetadeps" submodule.
    # Failure to do so is guaranteed to raise exceptions at BETSEE startup.
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    # For simplicity, py.test should remain the only hard dependency for
    # testing on local machines. While our setuptools-driven testing regime
    # optionally leverages third-party py.test plugins (e.g., "pytest-xdist"),
    # these plugins are *NOT* required for simple testing.
    #
    # A relatively modern version of py.test is required. Specifically:
    #
    # * At least version 3.7.0 or newer, which introduces the package scope for
    #   fixtures required to efficiently initialize and deinitialize
    #   application metadata singletons for unit tests. See also:
    #       https://docs.pytest.org/en/latest/fixture.html?highlight=scope#package-scope-experimental
    'pytest': '>= 3.7.0',
}
'''
Dictionary mapping from the :mod:`setuptools`-specific project name of each
mandatory testing dependency for this application to the suffix of a
:mod:`setuptools`-specific requirements string constraining this dependency.

See Also
----------
:data:`RUNTIME_MANDATORY`
    Further details on dictionary structure.
:download:`/doc/md/INSTALL.md`
    Human-readable list of these dependencies.
'''

# ....................{ LIBS ~ commands                   }....................
RequirementCommand = namedtuple('RequirementCommand', ('name', 'basename',))
RequirementCommand.__doc__ = '''
    Lightweight metadata describing a single external command required by an
    application dependency of arbitrary type (including optional, mandatory,
    runtime, testing, and otherwise).

    Attributes
    ----------
    name : str
        Human-readable name associated with this command (e.g., ``Graphviz``).
    basename : str
        Basename of this command to be searched for in the current ``${PATH}``.
    '''


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: Changes to this dictionary *MUST* be synchronized with:
# * Front-facing documentation (e.g., "doc/md/INSTALL.md").
# * Gitlab-CI configuration (e.g., the top-level "requirements-conda.txt" file).
# * Third-party platform-specific packages (e.g., Gentoo Linux ebuilds).
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REQUIREMENT_NAME_TO_COMMANDS = {
    'pydot': (RequirementCommand(name='Graphviz', basename='dot'),),
}
'''
Dictionary mapping from the :mod:`setuptools`-specific project name of each
application dependency (of any type, including optional, mandatory, runtime,
testing, or otherwise) requiring one or more external commands to a tuple of
:class:`RequirementCommand` instances describing these requirements.

See Also
----------
:download:`/doc/md/INSTALL.md`
    Human-readable list of these dependencies.
'''

# ....................{ GETTERS                           }....................
def get_runtime_mandatory_tuple() -> tuple:
    '''
    Tuple listing the :mod:`setuptools`-specific requirement string containing
    the mandatory name and optional version and extras constraints of each
    mandatory runtime dependency for this application, dynamically converted
    from the :data:`metadata.RUNTIME_MANDATORY` dictionary.
    '''

    # Avoid circular import dependencies.
    from betse.lib.setuptools import setuptool

    # Return this dictionary coerced into a tuple.
    return setuptool.get_requirements_str_from_dict(RUNTIME_MANDATORY)


def get_runtime_optional_tuple() -> tuple:
    '''
    Tuple listing the :mod:`setuptools`-specific requirement string containing
    the mandatory name and optional version and extras constraints of each
    optional runtime dependency for this application, dynamically converted
    from the :data:`metadata.RUNTIME_OPTIONAL` dictionary.
    '''

    # Avoid circular import dependencies.
    from betse.lib.setuptools import setuptool

    # Return this dictionary coerced into a tuple.
    return setuptool.get_requirements_str_from_dict(RUNTIME_OPTIONAL)


def get_testing_mandatory_tuple() -> tuple:
    '''
    Tuple listing the :mod:`setuptools`-specific requirement string containing
    the mandatory name and optional version and extras constraints of each
    mandatory testing dependency for this application, dynamically converted
    from the :data:`metadata.RUNTIME_OPTIONAL` dictionary.
    '''

    # Avoid circular import dependencies.
    from betse.lib.setuptools import setuptool

    # Return this dictionary coerced into a tuple.
    return setuptool.get_requirements_str_from_dict(TESTING_MANDATORY)
