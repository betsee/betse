#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2015 by Alexis Pietak & Cecil Curry
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

# ....................{ START                              }....................
# Import all constants defined by "betse.info" into the current namespace
# *BEFORE* subsequent logic possibly depending on the the version of the active
# Python interpreter, which such importation also validates.
#
# This awkward (albeit increasingly commonplace) snippet is required for
# reliable importation of metadata declared by and hence shared with the main
# codebase. Unfortunately, such metadata is *NOT* reliably importably via the
# conventional syntax (e.g., "from betse import info"). The reason is subtle.
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
with open('betse/info.py') as betse_info:
    exec(betse_info.read())

# ....................{ IMPORTS                            }....................
import setuptools

# ....................{ SETUP                              }....................
setuptools.setup(
    # ..................{ CORE                               }..................
    # Self-explanatory metadata. Since the "NAME" constant provided by
    # "betse.info" is uppercase *AND* since setuptools-installed package names
    # are commonly lowercase, such constant is coerced to lowercase.
    name = NAME.lower(),
    version = __version__,
    description = DESCRIPTION,
    author = AUTHORS,
    author_email='alexis.pietak@gmail.com',

    # PyPi-specific metadata.
    classifiers=[
        'Intended Audience :: Science/Research',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Topic :: Scientific/Engineering :: Visualization',
        # 'License :: ???',
    ],

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
