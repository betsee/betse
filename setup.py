#!/usr/bin/env python3
# Copyright 2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

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
    # Self-explanatory metadata.
    name = 'betse',
    version = '0.0.1',

    # List of all Python packages (i.e., directories containing zero or more
    # Python modules) to be installed. Currently, this includes the "betse"
    # package and all subpackages of such package excluding:
    #
    # * "betse.test" and all subpackages of such package, providing unit tests
    #   *NOT* intended to be installed with betse.
    packages = setuptools.find_packages(
        exclude=['betse.test', 'betse.test.*']),

    # Runtime dependencies.

    #FIXME: Formatted correctly? Verify, please.

    # Unit test-specific dependencies.
    tests_require = ['py.test'],
)

# --------------------( WASTELANDS                         )--------------------
