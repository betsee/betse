#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
:mod:`setuptools`-based makefile instrumenting all high-level administration
tasks (e.g., installation, freezing, test running) for this application.
'''

#FIXME: Replace this file (i.e., "setup.py") and the "requirements-pip.txt",
#"setup.cfg", "MANIFEST.in", and "Pipfile" files with the existing
#"pyproject.toml" file in concert with the third-party "poetry" project, which
#replaces Python's broken build management ecosystem (e.g., distutils,
#setuptools, pip, pipenv) with a single shell-friendly command patterned on
#industry-standard build management utilities published for other languages
#(e.g., Rust's "cargo"). We have manually inspected the "poetry" repository
#and, indeed, this is the setuptools killer we have long awaited.
#
#Note, however, that we may still require setuptools at runtime for dependency
#resolution -- which is quite alright, of course. Since everyone requires
#setuptools at installation-time, setuptools remains widely available and
#continuing to depend upon setuptools for a specific task remains sensible.

#FIXME: Consider wrapping on Windows with "pynsist", a framework for generating
#Windows installers bundling Python applications complete with a Python
#interpreter and requisite packages.

# ....................{ IMPORTS                           }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To avoid race conditions during setuptools-based installation, this
# module may import *ONLY* from packages guaranteed to exist at the start of
# installation. This includes all standard Python and application packages but
# *NOT* third-party dependencies, which if currently uninstalled will only be
# installed at some later time in the installation.
#
# Technically, this script may import from all subpackages and submodules of
# the this application's eponymous package. By Python mandate, the first
# element of the "sys.path" list is guaranteed to be the directory containing
# this script.  Python necessarily searches this directory for imports from the
# local version of this application *BEFORE* any other directories (including
# system directories containing older versions of this application). To quote:
#
#     "As initialized upon program startup, the first item of this list,
#      path[0], is the directory containing the script that was used to invoke
#      the Python interpreter."
#
# See also: https://stackoverflow.com/a/10097543/2809027
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

import setuptools
from betse import metadata
from betse.lib import libs
from betse.lib.setuptools.command import (
    supcmdfreeze, supcmdsymlink, supcmdtest)
from betse_setup import bupbuild, buputil

# ....................{ METADATA ~ seo                    }....................
_KEYWORDS = [
    'biology',
    'multiphysics',
    'science',
    'simulator',
    'disecrete exterior calculus',
]
'''
List of all lowercase alphabetic keywords synopsising this application.

These keywords may be arbitrarily selected so as to pretend to improve search
engine optimization (SEO). In actuality, they do absolutely nothing.
'''


# To minimize desynchronization woes, all
# "Programming Language :: Python :: "-prefixed strings are dynamically
# appended to this list by the init() function below.
_CLASSIFIERS = [
    # PyPI-specific version type. The number specified here is a magic constant
    # with no relation to this application's version numbering scheme. *sigh*
    'Development Status :: 5 - Production/Stable',

    # Sublist of all supported platform-specific CLI and GUI components.
    'Environment :: Console',
    'Environment :: MacOS X',
    'Environment :: Win32 (MS Windows)',
    'Environment :: X11 Applications',

    # Miscellaneous metadata.
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: BSD License',
    'Natural Language :: English',
    'Operating System :: OS Independent',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
]
'''
List of all PyPI-specific trove classifier strings synopsizing this
application.

Each such string *must* be contain either two or three `` :: `` substrings
delimiting human-readable capitalized English words formally recognized by the
:mod:`distutils`-specific ``register`` command.

See Also
----------
https://pypi.python.org/pypi?%3Aaction=list_classifiers
    Plaintext list of all trove classifier strings recognized by PyPI.
'''

# ....................{ OPTIONS                           }....................
# Setuptools-specific options. Keywords not explicitly recognized by either
# setuptools or distutils must be added to the above dictionary instead.
_setup_options = {
    # ..................{ CORE                              }..................
    # Self-explanatory metadata.
    'name':             metadata.PACKAGE_NAME,
    'version':          metadata.VERSION,
    'author':           metadata.AUTHORS,
    'author_email':     metadata.AUTHOR_EMAIL,
    'maintainer':       metadata.AUTHORS,
    'maintainer_email': metadata.AUTHOR_EMAIL,
    'description':      metadata.SYNOPSIS,
    'long_description': buputil.get_description(),
    'url':              metadata.URL_HOMEPAGE,
    'download_url':     metadata.URL_DOWNLOAD,

    # ..................{ PYPI                              }..................
    # PyPi-specific metadata.
    'classifiers': buputil.sanitize_classifiers(
        classifiers=_CLASSIFIERS,
        python_version_min_parts=metadata.PYTHON_VERSION_MIN_PARTS,
        python_version_minor_max=metadata.PYTHON_VERSION_MINOR_MAX,
    ),
    'keywords':    _KEYWORDS,
    'license':     metadata.LICENSE,

    # ..................{ DEPENDENCIES                      }..................
    # Mandatory runtime dependencies.
    'install_requires': libs.get_runtime_mandatory_tuple(),

    # Optional nuntime dependencies. Whereas mandatory dependencies are defined
    # as sequences, optional dependencies are defined as a dictionary mapping
    # from an arbitrary alphanumeric word to a sequence containing one or more
    # such dependencies. Such dependencies are then installable via "pip" by
    # suffixing the name of this project by the "["- and "]"-delimited key
    # defined below whose value lists the dependencies to be installed (e.g.,
    # "sudo pip3 install betse[all]", installing both the application and all
    # mandatory and optional dependencies required by the application).
    'extras_require': {
        # All optional runtime dependencies.
        'all': libs.get_runtime_optional_tuple(),
    },

    # Mandatory testing dependencies.
    'tests_require': libs.get_testing_mandatory_tuple(),

    # ..................{ PACKAGES                          }..................
    # List of the fully-qualified names of all Python packages (i.e.,
    # directories containing zero or more Python modules) to be installed,
    # including the top-level application package and all subpackages of that
    # package. This thus excludes:
    #
    # * The top-level test package and all subpackages of this package, test
    #   functionality *NOT* intended to be installed with this application.
    # * The top-level setup package and all subpackages of this package,
    #   setuptools functionality required only for application installation.
    # * "build", caching both setuptools metadata and a complete copy of this
    #   package, required only by a prior application installation.
    # * "freeze", providing PyInstaller-specific functionality required only for
    #   application freezing (i.e., conversion into an executable binary).
    #
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # WARNING: This inspection intentionally omits subdirectories containing no
    # "__init__.py" file, despite the remainder of the Python ecosystem
    # commonly accepting such subdirectories as subpackages.
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    'packages': setuptools.find_packages(exclude=(
        metadata.PACKAGE_NAME + '_test',
        metadata.PACKAGE_NAME + '_test.*',
        metadata.PACKAGE_NAME + '_setup',
        metadata.PACKAGE_NAME + '_setup.*',
        'build',
        'freeze',
    )),

    # ..................{ PATHS                             }..................
    # Cross-platform script wrappers dynamically created at installation time.
    'entry_points': {
        # CLI-specific scripts.
        'console_scripts': (
            '{} = {}.__main__:main'.format(
                metadata.SCRIPT_BASENAME, metadata.PACKAGE_NAME),
        ),
    },

    #FIXME; This isn't quite true. Undesirable files are excludable in this
    #file via the whitelist approach of "package_data" and/or blacklist approach
    #of "exclude_package_data". This functionality may or may not require the
    #external "setuptools-git" plugin, but is certainly feasible. See also the
    #comprehensive documentation at:
    #https://pythonhosted.org/setuptools/setuptools.html#including-data-files
    #FIXME: O.K.; if "'include_package_data': True" ends up not working for us,
    #contemplate:
    #
    #    'package_data': {
    #        # Include all files in all packages matching the following.
    #        '': ['*.txt', '*.xml', '*.special', '*.huh'],
    #    },

    # Install all data files (i.e., non-Python files) embedded in the Python
    # package tree for this application.
    #
    # Unlike Python packages, undesirable data files are excludable from
    # installation *ONLY* via the external "MANIFEST.in" file. This is
    # terrible, of course. (Did you expect otherwise?)
    #
    # Data files are *NOT* Python modules and hence should *NOT* be embedded in
    # the Python package tree. Sadly, the "data_files" key supported by
    # setuptools for this purpose is *NOT* cross-platform-portable and is thus
    # inherently broken. Why? Because this key either requires usage of
    # absolute paths *OR* relative paths relative to absolute paths defined by
    # "setup.cfg"; in either case, these paths are absolute. While the current
    # platform could be detected and the corresponding absolute path embedded
    # in 'data_files', that implementation would be inherently fragile. (That's
    # bad.) In lieu of sane setuptools support, we defer to the methodology
    # employed by everyone. Setuptools, your death is coming.
    'include_package_data': True,

    # Install to an uncompressed directory rather than a compressed archive.
    #
    # While nothing technically precludes the latter, doing so substantially
    # complicates runtime access of data files compressed into this archive
    # (e.g., with the pkg_resources.resource_filename() function). How so? By
    # decompressing this archive's contents into a temporary directory on
    # program startup and removing these contents on program shutdown. Since
    # there exists no guarantee this removal will actually be performed (e.g.,
    # due to preemptive SIGKILLs), compressed archives are inherently fragile.
    'zip_safe': False,
}
'''
Dictionary passed to the subsequent call to the :func:`setup` function.

This dictionary signifies the set of all application-specific :mod:`setuptools`
options. Submodules of the top-level :mod:`betse_setup` package subsequently
customize these options (e.g., by defining custom commands).
'''
# print('extras: {}'.format(setup_options['extras_require']))


# While currently empty, it's likely we'll want this again... someday.
_custom_metadata = {}
'''
Non-setuptools-specific metadata, used to inform custom subcommands (e.g.,
``freeze_file``) of supplementary metadata *not* already declared by the
:data:`setup_options` dictionary.

Setuptools raises fatal exceptions when the :data:`setup_options` dictionary
contains unrecognized keys. For safety, these keys are added to this dictionary
instead.
'''

# ....................{ COMMANDS                          }....................
# Define all application-specific setuptools commands.
for subcommand_submodule in bupbuild, supcmdfreeze, supcmdsymlink, supcmdtest:
    subcommand_submodule.add_subcommand(_setup_options, _custom_metadata)

# ....................{ SETUP                             }....................
setuptools.setup(**_setup_options)
