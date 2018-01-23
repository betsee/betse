#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
:mod:`setuptools`-based makefile instrumenting all high-level administration
tasks (e.g., installation, freezing, test running) for this application.
'''

#FIXME: Consider wrapping on Windows with "pynsist", a framework for generating
#Windows installers bundling Python applications complete with a Python
#interpreter and requisite packages.

# ....................{ IMPORTS                            }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To avoid race conditions during setuptools-based installation, this
# module may import *ONLY* from packages guaranteed to exist at the start of
# installation. This includes all standard Python and application packages but
# *NOT* third-party dependencies, which if currently uninstalled will only be
# installed at some later time in the installation.
#
# Technically, this script may import from all subpackages and submodules of the
# this application's eponymous package. By Python mandate, the first element of
# the "sys.path" list is guaranteed to be the directory containing this script.
# Python necessarily searches this directory for imports from the local version
# of this application *BEFORE* any other directories (including system
# directories containing older versions of this application). To quote:
#
#     "As initialized upon program startup, the first item of this list,
#      path[0], is the directory containing the script that was used to invoke
#      the Python interpreter."
#
# See also: https://stackoverflow.com/a/10097543/2809027
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

import setuptools
from betse import metadata
from betse.lib import libs
from betse_setup import build, freeze, symlink, test, buputil

# ....................{ METADATA                           }....................
# PyPI-specific metadata declared here rather than in the "betse.metadata"
# submodule, reducing space and time complexity during application startup. This
# metadata is relevant only to setuptools and hence irrelevant to the main
# codebase.

_PYTHON_VERSION_MINOR_MAX = 7
'''
Maximum minor stable version of this major version of Python currently released
(e.g., ``5`` if Python 3.5 is the most recent stable version of Python 3.x).
'''


_DESCRIPTION = None
'''
Human-readable multiline description of this application in reStructuredText
(reST) format.

To minimize synchronization woes, this description is identical to the contents
of the :doc:`/README.rst` file. When submitting this application package to
PyPI, this description is used verbatim as this package's front matter.
'''


# ....................{ METADATA ~ seo                     }....................
_KEYWORDS = ['biology', 'multiphysics', 'science', 'simulator',]
'''
List of all lowercase alphabetic keywords synopsising this application.

These keywords may be arbitrarily selected so as to pretend to improve search
engine optimization (SEO). In actuality, they do absolutely nothing.
'''


# To minimize desynchronization woes, all
# "Programming Language :: Python :: "-prefixed strings are dynamically appended
# to this list by the init() function below.
_CLASSIFIERS = [
    #FIXME: Replace with the following after the "1.0.0" release:
    #    'Development Status :: 5 - Production/Stable',

    # PyPI-specific version type. The number specified here is a magic constant
    # with no relation to this application's version numbering scheme. *sigh*
    'Development Status :: 4 - Beta',

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
List of all PyPI-specific trove classifier strings synopsizing this application.

Each such string *must* be contain either two or three `` :: `` substrings
delimiting human-readable capitalized English words formally recognized by the
:mod:`distutils`-specific ``register`` command.

See Also
----------
https://pypi.python.org/pypi?%3Aaction=list_classifiers
    Plaintext list of all trove classifier strings recognized by PyPI.
'''

# ....................{ INITIALIZERS                       }....................
def _init() -> None:
    '''
    Finalize the definition of all globals declared by this script.
    '''

    # Global variables assigned to below.
    global _DESCRIPTION

    # Major version of Python required by this application.
    PYTHON_VERSION_MAJOR = metadata.PYTHON_VERSION_MIN_PARTS[0]

    # Relative path of this application's front-facing documentation in
    # reStructuredText format, required by PyPI.
    DESCRIPTION_FILENAME = 'README.rst'

    # For each minor version of Python 3.x supported by this application,
    # formally classify this version as such.
    for python_version_minor in range(
        metadata.PYTHON_VERSION_MIN_PARTS[1], _PYTHON_VERSION_MINOR_MAX + 1):
        _CLASSIFIERS.append(
            'Programming Language :: Python :: {}.{}'.format(
                PYTHON_VERSION_MAJOR, python_version_minor,))
    # print('classifiers: {}'.format(_CLASSIFIERS))

    # Description read from this description file.
    try:
        _DESCRIPTION = buputil.get_chars(DESCRIPTION_FILENAME)
        # print('description: {}'.format(_DESCRIPTION))
    # If this file is *NOT* readable, print a non-fatal warning and reduce this
    # description to the empty string. While unfortunate, this description is
    # *NOT* required for most operations and hence mostly ignorable.
    except Exception as exception:
        _DESCRIPTION = ''
        buputil.output_warning(
            'Description file "{}" not found or not readable:\n{}'.format(
                DESCRIPTION_FILENAME, exception))

# ....................{ INITIALIZERS ~ main                }....................
# Finalize the definition of all globals declared by this module.
_init()

#FIXME: O.K.; it would appear that setuptools >= 38.0.0 now requires iterables
#to strictly be ordered sequences, whereas before it accepted unordered
#sequences. So, setuptools broke backwards compatibility yet again. Fortunately,
#the fix should be (mostly) trivial: just compact the
#"metadeps.RUNTIME_MANDATORY" dictionary passed to the "install_requires"
#parameter into a tuple of concatenated strings instead. We'll also need to do
#so for any other parameters accepting similar iterables. *sigh*

# print('mandatory runtime dependencies: {}'.format(metadeps.RUNTIME_MANDATORY))

# ....................{ OPTIONS                            }....................
# Setuptools-specific options. Keywords not explicitly recognized by either
# setuptools or distutils must be added to the above dictionary instead.
_setup_options = {
    # ..................{ CORE                               }..................
    # Self-explanatory metadata.
    'name':             metadata.PACKAGE_NAME,
    'version':          metadata.VERSION,
    'author':           metadata.AUTHORS,
    'author_email':     metadata.AUTHOR_EMAIL,
    'maintainer':       metadata.AUTHORS,
    'maintainer_email': metadata.AUTHOR_EMAIL,
    'description':      metadata.SYNOPSIS,
    'long_description': _DESCRIPTION,
    'url':              metadata.URL_HOMEPAGE,
    'download_url':     metadata.URL_DOWNLOAD,

    # ..................{ PYPI                               }..................
    # PyPi-specific metadata.
    'classifiers': _CLASSIFIERS,
    'keywords': _KEYWORDS,
    'license': metadata.LICENSE,

    # ..................{ DEPENDENCIES                       }..................
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
        'all': libs.get_runtime_optional_tuple(),
    },

    # Mandatory testing dependencies.
    'tests_require': libs.get_testing_mandatory_tuple(),

    # ..................{ PACKAGES                           }..................
    # List of all Python packages (i.e., directories containing zero or more
    # Python modules) to be installed. Currently, this includes the "betse"
    # package and all subpackages of this package excluding:
    #
    # * The top-level test package and all subpackages of this package, test
    #   functionality *NOT* intended to be installed with this application.
    # * The top-level setup package and all subpackages of this package,
    #   setuptools functionality required only for application installation.
    # * "build", caching both setuptools metadata and a complete copy of this
    #   package, required only by a prior application installation.
    # * "freeze", providing PyInstaller-specific functionality required only for
    #   application freezing (i.e., conversion into an executable binary).
    'packages': setuptools.find_packages(exclude=(
        metadata.PACKAGE_NAME + '_test',
        metadata.PACKAGE_NAME + '_test.*',
        metadata.PACKAGE_NAME + '_setup',
        metadata.PACKAGE_NAME + '_setup.*',
        'build',
        'freeze',
    )),

    # ..................{ PATHS                              }..................
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
    # inherently broken. Why? Because this key either requires usage of absolute
    # paths *OR* relative paths relative to absolute paths defined by
    # "setup.cfg"; in either case, these paths are absolute. While the current
    # platform could be detected and the corresponding absolute path embedded in
    # 'data_files', that implementation would be inherently fragile. (That's
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

    # ..................{ COMMANDS                           }..................
    # Set of all custom setuptools subcommands specific to this makefile (e.g.,
    # "sudo python3 setup.py symlink"), defaulting to the empty set. Each
    # subsequent call to the add_setup_commands() function iteratively performed
    # below adds one or more such subcommands to this set.
    'cmdclass': {},
}
'''
Dictionary passed to the subsequent call to the :func:`setup` function.

This dictionary signifies the set of all application-specific :mod:`setuptools`
options. Submodules of the top-level :mod:`betse_setup` package subsequently
customize these options (e.g., by defining custom commands).
'''
# print('extras: {}'.format(setup_options['extras_require']))


_setup_options_custom = {
    # While currently empty, it's likely we'll want this again... someday.
}
'''
Non-setuptools-specific metadata, used to inform custom subcommands (e.g.,
``freeze_file``) of other metadata *not* already declared by the
:data:`setup_options` dictionary.

Setuptools raises fatal exceptions when the :data:`setup_options` dictionary
contains unrecognized keys. For safety, these keys are added to this dictionary
instead.
'''

# ....................{ COMMANDS                           }....................
# Define all application-specific setuptools commands.
for setup_module in (build, freeze, symlink, test):
    setup_module.add_setup_commands(_setup_options_custom, _setup_options)

# ....................{ SETUP                              }....................
setuptools.setup(**_setup_options)
