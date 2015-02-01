#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''betse's setuptools-based makefile.'''

# ....................{ IMPORTS                            }....................
import setuptools

# ....................{ SETUP                              }....................
#FIXME; Add dependencies.
#FIXME; Add "pyside-uic" integration.
#FIXME; Add "py.test" integration. ("tox" as well, mayhaps?) For "py.test", see
#the following URL:
#    http://pytest.org/latest/goodpractises.html#integrating-with-distutils-python-setup-py-test
#Not terribly arduous, but we'll need to leverage some unctuous boilerplate.

setuptools.setup(
    # ..................{ CORE                               }..................
    # Self-explanatory metadata.
    name = 'betse',
    version = '0.0.1',
    description = (
        'betse (Bioelectric Tissue Simulation Environment) simulates '
        'propagation of electrical phenomena within biological tissue (e.g., '
        'ion channel-gated current flow).'),

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
