#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''`betse`'s `setuptools`-centric makefile.'''

# ....................{ IMPORTS                            }....................
import setuptools

# Import all constants defined by "betse.info" into the current namespace. This
# awkward (albeit increasingly commonplace) snippet is required for reliable
# importation of program metadata declared in the main codebase.  Unfortunately,
# such metadata is *NOT* reliably importably with conventional syntax (e.g., as
# "from betse import info"). The reasons, of course, are subtle.
#
# Importing packages in the main codebase implicitly imports the top-level
# "__init__.py" module provided by such codebase. If such module imports from at
# least one package *NOT* provided by stock Python installations (e.g., from
# packages installed as mandatory dependencies by this makefile), such
# importation will fail for users lacking such packages. While such module
# currently imports from no such packages, this race condition is sufficiently
# horrible as to warrant explicit circumvention: namely, by manually reading and
# evaluating the module defining such constants.
#
# This is horrible, but coding gets like that sometimes. We blame Guido.
with open('betse/info.py') as betse_info:
    exec(betse_info.read())

# ....................{ SETUP                              }....................
#FIXME; Add "pyside-uic" integration.

setuptools.setup(
    # ..................{ CORE                               }..................
    # Self-explanatory metadata.
    name = 'betse',
    version = __version__,
    description = DESCRIPTION,

    # ..................{ PATH                               }..................
    # List of all Python packages (i.e., directories containing zero or more
    # Python modules) to be installed. Currently, this includes the "betse"
    # package and all subpackages of such package excluding:
    #
    # * "betse.test" and all subpackages of such package, providing unit tests
    #   *NOT* intended to be installed with betse.
    # * "test", providing ad-hoc tests intended for developer use only.
    packages = setuptools.find_packages(
        exclude=['betse.test', 'betse.test.*', 'test', 'ui',]),

    #FIXME; This isn't quite true. Undesirable files are excludable via the
    #whitelist approach of "package_data" and/or blacklist approach of
    #"exclude_package_data". Such functionality may or may not require the
    #external "setuptools-git" plugin, but is certainly feasible. See also the
    #comprehensive documentation at:
    #https://pythonhosted.org/setuptools/setuptools.html#including-data-files

    # If True, *ALL* non-Python files will be installed as well. Unlike Python
    # packages, undesirable files are excludable from such installation via the
    # external "MANIFEST.in" file. (Inconvenient, but here we are.)
    include_package_data = True,

    # Install to an uncompressed directory rather than a compressed archive.
    # (While the latter could be made to work, we've yet to invest the effort.)
    zip_safe = False,

    # Cross-platform executable scripts dynamically created by setuptools at
    # installation time.
    entry_points = {
        # CLI-specific scripts.
        'console_scripts': ['betse = betse.cli.cli:main',],
        #FIXME: Create "betse.gui.gui".
        # GUI-specific scripts.
        'gui_scripts':  ['betse-qt = betse.gui.gui:main',],
    },

    # ..................{ DEPENDENCY                         }..................
    # Runtime dependencies. See "README.md".
    install_requires = [
        'numpy >= 1.9.0',
        'pyside >= 1.1.0',
        'pyyaml >= 3.10',
        'scipy >= 0.12.0',
        'matplotlib >= 1.3.0',
    ],

    # Unit test-specific dependencies. While such tests should also be runnable
    # under "py.test", "py.test" does *NOT* provide out-of-the-box support for
    # setuptools and hence is non-ideal.
    tests_require = ['nose >= 1.3.0'],

    # ..................{ TEST                               }..................
    # Name of the package running unit tests.
    test_suite = 'nose.collector',
)

# --------------------( WASTELANDS                         )--------------------
#FUXME; Add "py.test" integration. ("tox" as well, mayhaps?) For "py.test", see
#the following URL:
#    http://pytest.org/latest/goodpractises.html#integrating-with-distutils-python-setup-py-test
#Not terribly arduous, but we'll need to leverage some unctuous boilerplate.
