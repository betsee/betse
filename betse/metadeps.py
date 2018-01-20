#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Metadata constants synopsizing high-level application dependencies.
'''

# ....................{ IMPORTS                            }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To avoid race conditions during setuptools-based installation, this
# module may import *ONLY* from modules guaranteed to exist at the start of
# installation. This includes all standard Python and application modules but
# *NOT* third-party dependencies, which if currently uninstalled will only be
# installed at some later time in the installation.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

from betse.util.type import modules

# ....................{ LIBS ~ runtime : mandatory         }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: Changes to this subsection *MUST* be synchronized with:
# * Front-facing documentation (e.g., "doc/md/INSTALL.md").
# * The "betse.util.type.modules.DISTUTILS_PROJECT_NAME_TO_MODULE_NAME"
#   dictionary, converting between the setuptools-specific names listed below
#   and the Python-specific module names imported by this application.
# * Gitlab-CI configuration (e.g., the top-level "requirements-conda.txt" file).
# * Third-party platform-specific packages (e.g., Gentoo Linux ebuilds).
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RUNTIME_MANDATORY = {
    # setuptools is currently required at both install and runtime. At runtime,
    # setuptools is used to validate that dependencies are available.
    'setuptools': '>= 3.3',

    # Dependencies directly required by this application. Notably:
    #
    # * matplotlib >= 2.1.0 requires Numpy >= 1.7.1.
    # * SciPy >= 1.0.0 requires Numpy >= 1.8.2.
    #
    # SciPy wins.
    'Numpy':  '>= 1.8.2',
    'Pillow': '>= 2.3.0',
    'SciPy':  '>= 0.12.0',
    'dill':   '>= 0.2.3',

    # Matplotlib >= 1.5.0 is required for the newly added "viridis" colormap.
    'matplotlib': '>= 1.5.0',

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
  dependency, which may have no relation to the name of that project's top-level
  module or package (e.g., the ``PyYAML`` project's top-level package is
  :mod:`yaml`). For human readability in error messages, this name should
  typically be case-sensitively capitalized -- despite being parsed
  case-insensitively by :mod:`setuptools`.
* Value is either:
  * ``None`` or the empty string, in which case this dependency is unconstrained
    (i.e., any version of this dependency is sufficient).
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

# ....................{ LIBS ~ runtime : mandatory : yaml  }....................
RUNTIME_MANDATORY_YAML_PROJECT_NAME = None
'''
:mod:`setuptools`-specific project name of the current YAML implementation
satisfying this application's mandatory runtime YAML dependency.

This project is dynamically selected at submodule importation time based on the
most preferred implementation currently importable from the
:data:`RUNTIME_MANDATORY_YAML` dictionary.
'''


RUNTIME_MANDATORY_YAML = {
    'PyYAML': '>= 3.10',

    #FIXME: As "ruamel.yaml" now purports to safely roundtrip our documents,
    #let's give it a whirl and see how far we get plunge into the icy waters of
    #markup purgatory this time.

    # A relatively modern version of "ruamel.yaml" variants is required.
    # Specifically, this application requires:
    #
    # * At least version 0.15.24 or newer of "ruamel.yaml", which resolves a
    #   long-standing parser issue preventing overly complex YAML markup (such
    #   as ours) from being safely roundtripped:
    #   0.15.24 (2017-08-09):
    #   * (finally) fixed longstanding issue 23 (reported by Antony Sottile),
    #     now handling comment between block mapping key and value correctly
    # * The new "ruamel.yaml" API first introduced in 0.15.0. While older
    #   versions strictly adhere to the functional PyYAML-compatible API, newer
    #   versions break backward compatibility by entirely supplanting that API
    #   with a modern object-oriented approach. Supporting both isn't worth the
    #   substantial increase in maintenance debt.
    'ruamel.yaml': '>= 0.15.24',
}
'''
Dictionary mapping from the :mod:`setuptools`-specific project name of each
YAML implementation satisfying the same mandatory runtime dependency for this
application to the suffix of a :mod:`setuptools`-specific requirements string
constraining this dependency.

Since neither :mod:`distutils` nor :mod:`setuptools` support mutual exclusivity,
these alternatives are isolated into a separate dictionary implicitly merged
into the :data:`RUNTIME_MANDATORY` dictionary at submodule importation time.

Motivation
----------
Ideally, exactly one YAML implementation would be required. While
:mod:`ruamel.yaml` is the most preferred YAML implementation and hence the
obvious candidate, installation of that namespace package is sufficiently
non-trivial and frankly fragile to warrant falling back to less preferable but
more widely available alternatives.
'''

# ....................{ LIBS ~ runtime : optional          }....................
#FIXME: Should these be dependencies also be added to our "setup.py" metadata,
#perhaps as so-called "extras"? Contemplate. Consider. Devise.
RUNTIME_OPTIONAL = {
    # To simplify subsequent lookup at runtime, project names for optional
    # dependencies should be *STRICTLY LOWERCASE*. Since setuptools parses
    # project names case-insensitively, case is only of internal relevance.

    # Dependencies directly required by this application.
    'pympler': '>= 0.4.1',
    'ptpython': '>= 0.29',

    # A relatively modern version of NetworkX *EXCLUDING* 1.11, which
    # critically broke backwards compatibility by coercing use of the unofficial
    # inactive "pydotplus" PyDot fork rather than the official active "pydot"
    # PyDot project, is directly required by this application. NetworkX >= 1.12
    # reverted to supporting "pydot", thus warranting blacklisting of only
    # NetworkX 1.11. It is confusing, maybe?
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

# ....................{ LIBS ~ testing : mandatory         }....................
TESTING_MANDATORY = {
    # For simplicity, py.test should remain the only hard dependency for testing
    # on local machines. While our setuptools-driven testing regime optionally
    # leverages third-party py.test plugins (e.g., "pytest-xdist"), these
    # plugins are *NOT* required for simple testing.
    'pytest': '>= 2.5.0',
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

# ....................{ LIBS ~ command                     }....................
class DependencyCommand(object):
    '''
    Lightweight metadata describing a single external command required by some
    application dependency (of any type, including optional, mandatory, runtime,
    testing, or otherwise).

    Attributes
    ----------
    name : str
        Human-readable name associated with this command (e.g., ``Graphviz``).
    basename : str
        Basename of this command to be searched for in the current ``${PATH}``.
    '''

    def __init__(self, name: str, basename: str) -> None:
        self.name = name
        self.basename = basename


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: Changes to this dictionary *MUST* be synchronized with:
# * Front-facing documentation (e.g., "doc/md/INSTALL.md").
# * Gitlab-CI configuration (e.g., the top-level "requirements-conda.txt" file).
# * Third-party platform-specific packages (e.g., Gentoo Linux ebuilds).
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
EXTERNAL_COMMANDS = {
    'pydot': (DependencyCommand(name='Graphviz', basename='dot'),),
}
'''
Dictionary mapping from the :mod:`setuptools`-specific project name of each
application dependency (of any type, including optional, mandatory, runtime,
testing, or otherwise) requiring one or more external commands to a tuple of
:class:`DependencyCommand` instances describing these requirements.

See Also
----------
:download:`/doc/md/INSTALL.md`
    Human-readable list of these dependencies.
'''

# ....................{ INITIALIZERS                       }....................
def _init() -> None:
    '''
    Dynamically finalize the contents of the :data:`RUNTIME_MANDATORY` global
    declared by this submodule.

    Specifically, this function selects the first importable third-party YAML
    implementation (in descending order of preference) and merges the
    :mod:`setuptools`-specific metadata constraining this implementation from
    the :data:`RUNTIME_MANDATORY_YAML` global into the :data:`RUNTIME_MANDATORY`
    and :data:`RUNTIME_OPTIONAL` globals.

    Caveats
    ----------
    Note that this install- and runtime logic only applies to logic paths that
    actually execute this script. This excludes wheel installation, which avoids
    script execution for safety. (Because this is Python.)
    '''

    # Global variables assigned to below.
    global RUNTIME_MANDATORY_YAML_PROJECT_NAME

    # Setuptools-specific project name of the current YAML implementation
    # satisfying this application's mandatory runtime YAML dependency,
    # defaulting to the most preferred YAML framework. This has the beneficial
    # side effect of installing this framework at installation time if no other
    # YAML framework has been installed.
    RUNTIME_MANDATORY_YAML_PROJECT_NAME = 'ruamel.yaml'

    #FIXME: Handle version constraints as well. Since "setuptools" is itself a
    #non-standard dependency, handling such constraints in a portable manner
    #will probably mean leveraging "distutils" functionality for doing so.
    #Specifically, if a sufficiently new version of neither "ruamel.yaml" nor
    #"ruamel_yaml" is available, "yaml" should be fallen back to for safety.

    # Prefer "ruamel.yaml", the most actively maintained and hence preferred
    # YAML framework.
    if modules.is_module('ruamel.yaml'):
        pass
    # Fallback to "ruamel_yaml", the next most actively maintained and hence
    # preferred YAML framework. Unlike the more general-purpose "ruamel.yaml"
    # namespace package, "ruamel_yaml" is an Anaconda-specific non-namespace
    # package internally leveraged by the "conda" package manager but only
    # infrequently rebased against upstream changes.
    #FIXME: Uncomment once "ruamel_yaml" is sufficiently up-to-date.
    # elif modules.is_module('ruamel_yaml'):
    #     RUNTIME_MANDATORY_YAML_PROJECT_NAME = 'ruamel_yaml'
    # Fallback to PyYaml, the officially dead and hence least preferred YAML
    # framework.
    elif modules.is_module('yaml'):
        RUNTIME_MANDATORY_YAML_PROJECT_NAME = 'PyYAML'

    #FIXME: *DELETE THE FOLLOWING LINE AFTER* upstream "ruamel.yaml" roundtrip
    #issues are resolved. Currently, "ruamel.yaml" fails to properly roundtrip
    #our default simulation configuration and hence must *NOT* be used. *sigh*
    #
    #After upstream both resolves these issues *AND* releases a new stable
    #release, modify the "ruamel.yaml" version constraints above to reflect the
    #new minimum required version. Naturally, this implies we won't be enabling
    #"ruamel.yaml" support anytime soon.
    #FIXME: *UGH.* "ruamel.yaml" appears to have a core deficiency with respect
    #to map key comments: specifically, the roundtripper silently produces
    #malformed YAML when encountering map key comments. This is quite obviously
    #unacceptable for our heavily commented YAML files. Sadly, upstream has
    #known about this open issue for several years -- implying a fix is distant
    #at beast. Until resolved, "ruamel.yaml" *MUST* be ignored. For discussion,
    #see the following upstream issue:
    #    https://bitbucket.org/ruamel/yaml/issues/146/torture-test-breaks-the-fragile-back-of
    RUNTIME_MANDATORY_YAML_PROJECT_NAME = 'PyYAML'

    # Enforce installation of the preferred YAML framework detected above.
    RUNTIME_MANDATORY[RUNTIME_MANDATORY_YAML_PROJECT_NAME] = (
        RUNTIME_MANDATORY_YAML[RUNTIME_MANDATORY_YAML_PROJECT_NAME])

    # Merge the requirement strings for all YAML frameworks into those for all
    # optional dependencies, permitting the former to be treated like the latter
    # (e.g., by the betse.lib.libs.die_unless_runtime_optional() function).
    RUNTIME_OPTIONAL.update(RUNTIME_MANDATORY_YAML)

# ....................{ MAIN                               }....................
# Initialize this submodule.
_init()
