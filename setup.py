#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
BETSE's `setuptools`-based makefile.
'''

#FIXME; Add "pyside-uic" integration. This is feasible as demonstrated by the
#following URL, which appears to be the only online reference to such practice.
#We could leverage such logic by defining a new "setup_pyside.py" file in the
#same directory as this file containing the class defined by:
#
#   https://gist.github.com/ivanalejandro0/6758741
#
#We'll want to retain the "pyside-rcc"-specific portion of that code as well,
#which compiles ".qrc" resource files (...in XML format, of course) to ".py"
#files. See the useful answer here concerning usage from the Python perspective:
#
#    https://stackoverflow.com/questions/22508491/a-py-file-which-compiled-from-qrc-file-using-pyside-rcc-does-not-work

#FIXME: Specify all setup() metadata keys listed here:
#https://docs.python.org/3/distutils/setupscript.html#additional-meta-data

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
from betse_setup import build, freeze, symlink, test

# ....................{ OPTIONS                            }....................
# Setuptools-specific options. Keywords not explicitly recognized by either
# setuptools or distutils must be added to the above dictionary instead.
setup_options = {
    # ..................{ CORE                               }..................
    # Self-explanatory metadata. Since the "NAME" constant provided by
    # "betse.info" is uppercase *AND* since setuptools-installed package names
    # are commonly lowercase, such constant is coerced to lowercase.
    'name': metadata.NAME.lower(),
    'version': metadata.__version__,
    'description': metadata.DESCRIPTION,
    'author': metadata.AUTHORS,
    'author_email': 'alexis.pietak@gmail.com',
    'url': 'https://gitlab.com/betse/betse',

    # PyPi-specific metadata.
    'keywords': 'science research visualization',
    'classifiers': [
        'Intended Audience :: Science/Research',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Topic :: Scientific/Engineering :: Visualization',
        # 'License :: ???',
    ],

    # ..................{ COMMAND                            }..................
    # Set of all custom setuptools subcommands specific to this makefile (e.g.,
    # "sudo python3 setup.py symlink"), defaulting to the empty set. Each
    # subsequent call to the add_setup_commands() function iteratively performed
    # below adds one or more such subcommands to this set.
    'cmdclass': {},

    # ..................{ PATH                               }..................
    # Cross-platform script wrappers dynamically created at installation time.
    'entry_points': {
        # CLI-specific scripts.
        'console_scripts': [
            metadata.SCRIPT_NAME_CLI + ' = betse.cli.__main__:main',
        ],
        # 'console_scripts': [SCRIPT_NAME_CLI + ' = betse.cli.clicli:main',],

        #FIXME: After creating a BETSE GUI, uncomment the following logic.
        # GUI-specific scripts.
        #'gui_scripts':  [
        #     metadata.SCRIPT_NAME_GUI + ' = betse.gui.guicli:main',
        #],
    },

    # List of all Python packages (i.e., directories containing zero or more
    # Python modules) to be installed. Currently, this includes the "betse"
    # package and all subpackages of this package excluding:
    #
    # * "betse_test" and all subpackages of this package, providing
    #   test-specific functionality *NOT* intended to be installed with BETSE.
    # * "betse_setup" and all subpackages of this package, providing
    #   setuptools-specific functionality required only for BETSE installation.
    # * "freeze", providing PyInstaller-specific functionality required only for
    #   BETSE freezing (i.e., conversion into an executable binary).
    # * "old", providing obsolete code to be shortly put to pasture.
    'packages': setuptools.find_packages(
        exclude = [
            # 'betse.data', 'betse.data.*',
            'betse_test', 'betse_test.*',
            'betse_setup', 'betse_setup.*',
            'freeze',
            'old',
        ],
    ),

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
    # the Python package tree. Sadly, the 'data_files' key supported by
    # setuptools for this purpose is *NOT* cross-platform-portable and hence
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

    # ..................{ DEPENDENCY                         }..................
    # Runtime dependencies.
    'install_requires': metadata.DEPENDENCIES_RUNTIME,

    #FIXME; Switch to "py.test". Setuptools integration isn't terribly arduous,
    #but does require some unctuous boilerplate:
    #    https://pytest.org/latest/goodpractices.html

    # Testing dependencies.
    'tests_require': metadata.DEPENDENCIES_TESTING,

    # ..................{ TEST                               }..................
    # Name of the package running unit tests.
    # 'test_suite': 'nose.collector',
}
'''
Dictionary passed to the subsequent call to `setup()`.

This dictionary signifies the set of all `betse`-specific `setuptools` options.
Modules in the `betse`-specific `setup` package customize such options (e.g., by
defining custom commands).
'''


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
