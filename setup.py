#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''`betse`'s `setuptools`-based makefile.'''

#FIXME: Define a new "symlink" setuptools command, strongly inspired by the
#IPython command of the same name. See:
#    https://github.com/ipython/ipython/blob/master/setupbase.py
#    https://github.com/ipython/ipython/blob/master/setup.py
#To support such cruft, we probably want a new setuptools-specific package tree
#(e.g., a new top-level directory "setup" containing at least files
#"__init__.py" and "symlink.py", the latter implementing the "symlink" and
#"unsymlink" commands).

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

# ....................{ START                              }....................
# Import all constants defined by "betse.metadata" into the current namespace
# *BEFORE* subsequent logic possibly depending on the the version of the active
# Python interpreter, which such importation also validates.
#
# This awkward (albeit increasingly commonplace) snippet is required for
# reliable importation of metadata declared by and hence shared with the main
# codebase. Unfortunately, such metadata is *NOT* reliably importably via the
# conventional syntax (e.g., "from betse import metadata"). The reason is
# subtle.
#
# Importing packages from the main codebase implicitly imports such codebase's
# top-level "__init__.py" module. If such module imports from at least one
# package *NOT* provided by stock Python installations (e.g., from packages
# installed as mandatory dependencies by this makefile), such importation will
# fail for users lacking such packages. While such module currently imports from
# no such packages, this race condition is sufficiently horrible as to warrant
# explicit circumvention: namely, by manually reading and evaluating the module
# defining such constants.
#
# This is horrible, but coding gets like that sometimes. We blame Guido.
with open('betse/metadata.py') as betse_metadata:
    exec(betse_metadata.read())

# ....................{ IMPORTS                            }....................
from setup import scripts, symlink
import setuptools

# ....................{ OPTIONS                            }....................
setup_options = {
    # ..................{ CORE                               }..................
    # Self-explanatory metadata. Since the "NAME" constant provided by
    # "betse.info" is uppercase *AND* since setuptools-installed package names
    # are commonly lowercase, such constant is coerced to lowercase.
    'name': NAME.lower(),
    'version': __version__,
    'description': DESCRIPTION,
    'author': AUTHORS,
    'author_email': 'alexis.pietak@gmail.com',

    # PyPi-specific metadata.
    'classifiers': [
        'Intended Audience :: Science/Research',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Topic :: Scientific/Engineering :: Visualization',
        # 'License :: ???',
    ],

    # ..................{ COMMAND                            }..................
    # Custom commands specific to this makefile called in the customary way
    # (e.g., "sudo python3 setup.py symlink").
    'cmdclass': {},

    # ..................{ PATH                               }..................
    # List of all Python packages (i.e., directories containing zero or more
    # Python modules) to be installed. Currently, this includes the "betse"
    # package and all subpackages of such package excluding:
    #
    # * "betse.test" and all subpackages of such package, providing unit tests
    #   *NOT* intended to be installed with betse.
    # * "setup", providing setuptools-specific packages.
    # * "test", providing ad-hoc tests intended for developer use only.
    'packages': setuptools.find_packages(
        exclude = ['betse.test', 'betse.test.*', 'setup', 'test', 'ui',],
    ),

    #FIXME; This isn't quite true. Undesirable files are excludable via the
    #whitelist approach of "package_data" and/or blacklist approach of
    #"exclude_package_data". Such functionality may or may not require the
    #external "setuptools-git" plugin, but is certainly feasible. See also the
    #comprehensive documentation at:
    #https://pythonhosted.org/setuptools/setuptools.html#including-data-files

    # If True, *ALL* non-Python files will be installed as well. Unlike Python
    # packages, undesirable files are excludable from such installation via the
    # external "MANIFEST.in" file. (Inconvenient, but here we are.)
    'include_package_data': True,

    # Install to an uncompressed directory rather than a compressed archive.
    # (While the latter could be made to work, we've yet to invest the effort.)
    'zip_safe': False,

    # Cross-platform executable scripts dynamically created by setuptools at
    # both installation and symlink time.
    'entry_points': {
        # CLI-specific scripts.
        'console_scripts': [SCRIPT_NAME_CLI + ' = betse.cli.clicli:main',],

        #FIXME: Create "betse.gui.guicli".
        # GUI-specific scripts.
        'gui_scripts':  [SCRIPT_NAME_GUI + ' = betse.gui.guicli:main',],
    },

    # ..................{ DEPENDENCY                         }..................
    # Runtime dependencies. See "README.md".
    'install_requires': REQUIREMENTS,

    # Unit test-specific dependencies. While such tests should also be runnable
    # under "py.test", "py.test" does *NOT* provide out-of-the-box support for
    # setuptools and hence is non-ideal.
    'tests_require': ['nose >= 1.3.0'],

    # ..................{ TEST                               }..................
    # Name of the package running unit tests.
    'test_suite': 'nose.collector',
}
'''
Dictionary passed to the subsequent call to `setup()`.

This dictionary signifies the set of all `betse`-specific `setuptools` options.
Modules in the `betse`-specific `setup` package customize such options (e.g., by
defining custom commands).
'''

# ....................{ COMMANDS                           }....................
# Define all BETSE-specific setuptools commands.
for setup_module in (scripts, symlink):
    setup_module.add_commands(setup_options)

# ....................{ SETUP                              }....................
setuptools.setup(**setup_options)

# --------------------( WASTELANDS                         )--------------------
    # 'install_requires': [
    #     'numpy >= 1.9.0',
    #     'pyside >= 1.1.0',
    #     'pyyaml >= 3.10',
    #     'scipy >= 0.12.0',
    #     'matplotlib >= 1.3.0',
    # ],

# ENTRY_POINTS = {
#     # CLI-specific scripts.
#     'console_scripts': ['betse = betse.cli.cli:main',],
#
#     #FUXME: Create "betse.gui.gui".
#     # GUI-specific scripts.
#     'gui_scripts':  ['betse-qt = betse.gui.gui:main',],
# }
# '''
# Dictionary defining cross-platform executable scripts dynamically created by
# `setuptools` at both installation and symlink time.
#
# For `setuptools` compatibility, this dictionary assumes the same structure as
# the `entry_points` field passed to `setup()`.
# '''

#FUXME; Add "py.test" integration. ("tox" as well, mayhaps?) For "py.test", see
#the following URL:
#    http://pytest.org/latest/goodpractises.html#integrating-with-distutils-python-setup-py-test
#Not terribly arduous, but we'll need to leverage some unctuous boilerplate.
