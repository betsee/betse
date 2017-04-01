#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
BETSE's `setuptools`-based makefile.
'''

#FIXME; Add "pyside-uic" integration. This is feasible as demonstrated by the
#following URL, which appears to be the only online reference for this practice.
#We could leverage this logic by defining a new "setup_pyside.py" file in the
#same directory as this file containing the class defined by:
#
#   https://gist.github.com/ivanalejandro0/6758741
#
#We'll want to retain the "pyside-rcc"-specific portion of that code as well,
#which compiles ".qrc" resource files (...in XML format, of course) to ".py"
#files. See the useful answer here concerning usage from the Python perspective:
#
#    https://stackoverflow.com/questions/22508491/a-py-file-which-compiled-from-qrc-file-using-pyside-rcc-does-not-work

# ....................{ IMPORTS                            }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To avoid race conditions during setuptools-based installation, this
# module may import *ONLY* from packages guaranteed to exist at the start of
# installation. This includes all standard Python and BETSE packages but *NOT*
# third-party dependencies, which if currently uninstalled will only be
# installed at some later time in the installation.
#
# Technically, this script may import from all packages in the BETSE codebase
# including the top-level "betse", "betse_setup", and "betse_test" packages.
# By Python mandate, the first element of "sys.path" is guaranteed to be the
# directory containing this script. Hence, Python necessarily searches this
# directory for imports from the local version of BETSE *BEFORE* any other
# directories (including system directories containing previously installed
# versions of BETSE). To quote:
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
from betse.util.io import stderrs
from betse.util.path import files
from betse_setup import build, freeze, symlink, test

# ....................{ METADATA                           }....................
# PyPI-specific metadata declared here rather than in "betse.metadata" to reduce
# both space and time complexity for BETSE startup. This metadata is effectively
# setuptools-specific and hence irrelevant to the main codebase.

_PYTHON_VERSION_MINOR_MAX = 6
'''
Maximum minor stable version of this major version of Python currently released
(e.g., `5` if Python 3.5 is the most recent stable version of Python 3.x).
'''


_DESCRIPTION = None
'''
Human-readable multiline description of this application in reStructuredText
(reST) format.

To minimize synchronization woes, this description is identical to the contents
of the `README.rst` file. When submitting this application package to PyPI,
this description is used verbatim as this package's front matter.
'''


_KEYWORDS = ['biology', 'multiphysics', 'science', 'simulator',]
'''
List of all lowercase alphabetic keywords synopsising this application.

These keywords may be arbitrarily selected so as to pretend to improve search
engine optimization (SEO). In actuality, they do absolutely nothing.
'''

# ....................{ METADATA ~ trove                   }....................
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

Each such string _must_ be contain either two or three ` :: ` substrings
delimiting human-readable capitalized English words formally recognized by the
`distutils`-specific `register` command.

See Also
----------
https://pypi.python.org/pypi?%3Aaction=list_classifiers
    Plaintext list of all trove classifier strings recognized by PyPI.
'''

# ....................{ INITIALIZERS                       }....................
def init() -> None:
    '''
    Finalize the definition of all globals declared by this module.
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
        metadata.PYTHON_VERSION_MIN_PARTS[1],
        _PYTHON_VERSION_MINOR_MAX + 1):
        _CLASSIFIERS.append(
            'Programming Language :: Python :: {}.{}'.format(
                PYTHON_VERSION_MAJOR, python_version_minor,))
    # print('classifiers: {}'.format(_CLASSIFIERS))

    # Description read from this description file.
    try:
        _DESCRIPTION = files.get_chars(DESCRIPTION_FILENAME)
    # If this file is *NOT* readable, print a non-fatal warning and reduce this
    # description to the empty string. While unfortunate, this description is
    # *NOT* required for most operations and hence mostly ignorable.
    except Exception as exception:
        _DESCRIPTION = ''
        stderrs.output(
            'Description file "{}" not found or not readable:\n{}'.format(
                DESCRIPTION_FILENAME, exception))


# Finalize the definition of all globals declared by this module.
init()

# ....................{ OPTIONS                            }....................
# Setuptools-specific options. Keywords not explicitly recognized by either
# setuptools or distutils must be added to the above dictionary instead.
setup_options = {
    # ..................{ CORE                               }..................
    # Self-explanatory metadata.
    'name':             metadata.PACKAGE_NAME,
    'version':          metadata.__version__,
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
    # Mandatory nuntime dependencies.
    'install_requires': metadata.DEPENDENCIES_RUNTIME_MANDATORY,

    # Optional nuntime dependencies. Whereas mandatory dependencies are defined
    # as sequences, optional dependencies are defined as a dictionary mapping
    # from an arbitrary alphanumeric word to a sequence containing one or more
    # such dependencies. Such dependencies are then installable via "pip" by
    # suffixing the name of this project by the "["- and "]"-delimited key
    # defined below whose value lists the dependencies to be installed (e.g.,
    # "sudo pip3 install betse[all]", installing both BETSE and all mandatory
    # and optional dependencies transitively required by BETSE).
    'extras_require': {
        # All optional runtime dependencies. Since the
        # "DEPENDENCIES_RUNTIME_OPTIONAL" global is a dictionary rather than a
        # sequence, a function converting this global into a tuple is called.
        'all': libs.get_runtime_optional_tuple(),
    },

    # Mandatory testing dependencies.
    'tests_require': metadata.DEPENDENCIES_TESTING_MANDATORY,

    # ..................{ PACKAGES                           }..................
    # List of all Python packages (i.e., directories containing zero or more
    # Python modules) to be installed. Currently, this includes the "betse"
    # package and all subpackages of this package excluding:
    #
    # * "betse_test" and all subpackages of this package, providing
    #   test-specific functionality *NOT* intended to be installed with BETSE.
    # * "betse_setup" and all subpackages of this package, providing
    #   setuptools-specific functionality required only for BETSE installation.
    # * "build", caching both setuptools-specific metadata and a complete copy
    #   of this package, required only by a prior BETSE installation.
    # * "freeze", providing PyInstaller-specific functionality required only for
    #   BETSE freezing (i.e., conversion into an executable binary).
    'packages': setuptools.find_packages(
        exclude = [
            'betse_test', 'betse_test.*',
            'betse_setup', 'betse_setup.*',
            'build',
            'freeze',
        ],
    ),

    # ..................{ PATHS                              }..................
    # Cross-platform script wrappers dynamically created at installation time.
    'entry_points': {
        # CLI-specific scripts.
        'console_scripts': [
            metadata.SCRIPT_NAME_CLI + ' = betse.__main__:main',
        ],

        #FIXME: After creating a BETSE GUI, uncomment the following logic.
        # GUI-specific scripts.
        #'gui_scripts':  [
        #     metadata.SCRIPT_NAME_GUI + ' = betse.gui.guicli:main',
        #],
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
Dictionary passed to the subsequent call to `setup()`.

This dictionary signifies the set of all `betse`-specific `setuptools` options.
Modules in the `betse`-specific `setup` package customize such options (e.g., by
defining custom commands).
'''
# print('extras: {}'.format(setup_options['extras_require']))


setup_options_custom = {
    # While currently empty, it's likely we'll want this again... someday.
}
'''
Non-setuptools-specific metadata, used to inform custom subcommands (e.g.,
`freeze_file`) of other metadata _not_ already declared by the `setup_options`
dictionary.

Setuptools raises fatal exceptions when the `setup_options` dictionary contains
unrecognized keys. For safety, these keys are added to this dictionary instead.
'''

# ....................{ COMMANDS                           }....................
# Define all BETSE-specific setuptools commands.
for setup_module in (build, freeze, symlink, test):
    setup_module.add_setup_commands(setup_options_custom, setup_options)

# ....................{ SETUP                              }....................
setuptools.setup(**setup_options)
