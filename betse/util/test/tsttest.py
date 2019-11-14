#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level **test suite** (e.g., collection of one or more functional or unit
tests collectively exercising problematic features of this application)
facilities.
'''

# ....................{ IMPORTS                           }....................
# from betse.util.io.log import logs
# from betse.util.type.types import type_check
from betse.util.type.decorator.decmemo import func_cached

# ....................{ TESTERS                           }....................
@func_cached
def is_testing() -> bool:
    '''
    ``True`` only if the active Python interpreter is currently exercising
    tests (e.g., via the :mod:`pytest` test harness).

    Caveats
    ----------
    **This function is only safely callable after application initialization,**
    implying that this function should *never* be called from module scope.
    For generality, this function requires the application metadata singleton
    to have already been defined.

    Raises
    ----------
    BetseMetaAppException
        If this application is uninitialized (i.e., the application metadata
        singleton is currently undefined).
    '''

    # There are numerous means of implementing this test, including:
    #
    # * Testing whether or not the "pytest" package has been previously
    #   imported. This is the most commonly recommended approach on both
    #   StackOverflow and in Django commentary. While trivial to implement,
    #   this approach has the distinct disadvantage of generating spurious
    #   false positives. Avoiding such issues requires *NEVER* importing from
    #   the "pytest" package from the main codebase. To generalize testing
    #   functionality across this application and downstream consumers (e.g.,
    #   BETSEE), the "betse.util.test" subpackage commonly does just that.
    #   Ergo, this approach is effectively infeasible here.
    # * Enabling a private boolean constant (e.g.,
    #   "betse.util.test.tests._IS_TESTING") at py.test session start time
    #   (e.g., from a session-wide autouse fixture or
    #   betse_test.conftest.pytest_configure() hook). While technically
    #   feasible, this approach is considerably less trivial *AND* invites
    #   subtle complications at chicken-and-egg application initialization
    #   time. Quite simply: the rewards are low and the prices are high.
    #
    # The current approach circumvents all of the above issues. Like the first
    # approach above, this approach is trivial to implement; unlike the first
    # approach above, this approach does *NOT* generate spurious false
    # positives. Why? Because the main codebase is guaranteed to *NEVER* import
    # from this application's test suite. Doing so would fundamentally violate
    # sanity in numerous ways and, in any case, is never desirable. In
    # particular, this application's test suite is *NOT* bundled with this
    # application and hence *ONLY* importable during development.

    # Avoid circular import dependencies.
    from betse.util.app.meta import appmetaone
    from betse.util.py.module import pymodname

    # Name of the root package of this application's test suite.
    test_package_name = appmetaone.get_app_meta().test_package_name

    # Return true only if that package has been previously imported from.
    return pymodname.is_imported(test_package_name)


#FIXME: The ideal implementation would detect whether the basename of the
#command invoking the parent process is "tox" or not. While doing so under
#Linux specifically is trivial, generalizing this detection to both macOS and
#Windows is highly non-trivial; doing so effectively requires the third-party
#"psutil" dependency, which is non-ideal. Nonetheless, perhaps it is time to
#bite that bullet and simply do so. We would like to leverage "psutil"
#elsewhere for various purposes (e.g., detecting PowerShell). Or... perhaps
#not. Consider this heuristic:
#
#* If "psutil" is importable, leverage the solution given at:
#  https://stackoverflow.com/a/24115041/2809027
#  This is guaranteed to be the most portable and reliable approach. Note,
#  however, that there's probably a slightly simpler implementation: the
#  psutil.Process.parent() function, obviating the need for PIDs. (Noice!) See:
#  https://github.com/giampaolo/psutil#process-management
#* Else:
#  * Under Linux, leverage the solution given at:
#    https://stackoverflow.com/a/24114907/2809027
#  * Under Windows, leverage the solution already documented at the
#    "betse.util.os.brand.windows" submodule. For now, we could simply raise an
#    exception.
#  * Under macOS, raise an exception. Ain't no one got time for that.
@func_cached
def is_tox() -> bool:
    '''
    ``True`` only if the active Python interpreter is currently exercising
    tests via :mod:`tox`, a high-level utility isolating lower-level test
    harnesses (e.g., :mod:`nose`, :mod:`pytest`) to test-specific cached venvs.
    '''

    # Avoid circular import dependencies.
    # from betse.util.py import pyvenv

    #FIXME: Ugh! Ain't got no idea how to actually detect "tox" reliably.
    #Probably have to check the process list for a parent named "tox". *sigh*
    return False

    # Return true only if...
    # return (
    #     # The active Python interpreter is isolated to a venv, which "tox"
    #     # guarantees, *AND*...
    #     pyvenv.is_venv() and
    # )
