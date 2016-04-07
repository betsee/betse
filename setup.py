#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''`betse`'s `setuptools`-based makefile.'''

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
#https://docs.python.org/2/distutils/setupscript.html#additional-meta-data

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
from setup import build, freeze, symlink
import setuptools

# ....................{ OPTIONS                            }....................
# Non-setuptools-specific metadata, used to inform custom subcommands (e.g.,
# "freeze_file") of other metadata *NOT* already declared by the "setup_options"
# dictionary defined below. Since setuptools raises fatal exceptions on such
# dictionary containing unrecognized keys, such keys are collected here.
metadata = {
    # While currently empty, it's likely we'll want this again... someday.
}

# Setuptools-specific options. Keywords not explicitly recognized by either
# setuptools or distutils must be added to the above dictionary instead.
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
    'url': 'http://www.mindshines.com',

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
    # Custom commands specific to this makefile called in the customary way
    # (e.g., "sudo python3 setup.py symlink").
    'cmdclass': {},

    # ..................{ PATH                               }..................
    # Cross-platform script wrappers dynamically created at installation time.
    'entry_points': {
        # CLI-specific scripts.
        'console_scripts': [SCRIPT_NAME_CLI + ' = betse.cli.__main__:main',],
        # 'console_scripts': [SCRIPT_NAME_CLI + ' = betse.cli.clicli:main',],

        #FIXME: Create "betse.gui.guicli".
        # GUI-specific scripts.
        'gui_scripts':  [SCRIPT_NAME_GUI + ' = betse.gui.guicli:main',],
    },

    # List of all Python packages (i.e., directories containing zero or more
    # Python modules) to be installed. Currently, this includes the "betse"
    # package and all subpackages of such package excluding:
    #
    # * "betse.test" and all subpackages of such package, providing unit tests
    #   *NOT* intended to be installed with betse.
    # * "setup", providing setuptools-specific packages.
    # * "test", providing ad-hoc tests intended for developer use only.
    'packages': setuptools.find_packages(
        exclude = [
            # 'betse.data', 'betse.data.*',
            'betse.test', 'betse.test.*',
            'setup',
            'test',
            'ui',
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

    # Install to uncompressed directories rather than compressed archives.
    #
    # While nothing technically precludes the latter, doing so substantially
    # complicates runtime access of data files compressed into such archives.
    # Specifically, doing so complicates usage of the
    # pkg_resources.resource_filename() function by requiring runtime
    # decompression of such archive's contents to a temporary directory and
    # removal of such contents at program exit. There is no guarantee that such
    # removal will actually be run (e.g., due to preemptive SIGKILLs), implying
    # such approach to be inherently fragile and hardly worth the effort.
    'zip_safe': False,

    # ..................{ DEPENDENCY                         }..................
    # Runtime dependencies.
    'install_requires': DEPENDENCIES_RUNTIME,

    #FIXME; Switch to "py.test". Setuptools integration isn't terribly arduous,
    #but does require some unctuous boilerplate:
    #    https://pytest.org/latest/goodpractices.html

    # Testing dependencies.
    'tests_require': DEPENDENCIES_TESTING,

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

# ....................{ COMMANDS                           }....................
# Define all BETSE-specific setuptools commands.
for setup_module in (build, freeze, symlink):
    setup_module.add_setup_commands(metadata, setup_options)

# ....................{ SETUP                              }....................
setuptools.setup(**setup_options)
