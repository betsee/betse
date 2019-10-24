#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Fixtures and fixture classes creating temporary simulation configurations
isolated to specific tests, which typically modify the contents of these
configurations so as to exercise specific feature sets and edge cases.
'''

# ....................{ IMPORTS                           }....................
from betse_test.fixture import initter
from betse_test.fixture.simconf.simconfclser import SimConfTestInternal
from pytest import fixture
from py._path.local import LocalPath

# ....................{ FIXTURES                          }....................
@fixture
def betse_sim_conf_default(betse_temp_dir: LocalPath) -> SimConfTestInternal:
    '''
    Per-test fixture creating a temporary default simulation configuration file
    and returning a wrapper around this file.

    Unlike the minified simulation configuration created by the
    :func:`betse_sim_conf` fixture and leveraged by most tests, the default
    simulation configuration created by this fixture remains unmodified (except
    for unavoidably disabling interactive simulation features, which
    non-interactive testing requires). Tests leveraging this fixture incur a
    significant performance penalty but can expose edge-case issues obscured by
    minification, including computational instability produced by the default
    non-minified time steps.

    Parameters
    ----------
    betse_temp_dir : LocalPath
        Wrapper around a temporary directory isolated to the current test.

    Returns
    ----------
    SimConfTestInternal
        Wrapper around a temporary default simulation configuration file.

    See Also
    ----------
    :func:`betse_sim_conf`
        Further details, ignoring minification performed by this fixture.
    '''

    # Defer heavyweight imports.
    from betse.science.parameters import Parameters

    # Initialize the application metadata singleton, as required by the
    # subsequent access of the "Parameters.conf_default_filename" class
    # property.
    initter.init_app()

    # Wrapper wrapping the default simulation configuration file copied into
    # this temporary directory and sanitized therein.
    sim_state = SimConfTestInternal(
        src_conf_filename=Parameters.conf_default_filename,
        trg_conf_filepath=betse_temp_dir.join('sim_config.yaml'))

    # Return this wrapper *WITHOUT* minifying this configuration.
    return sim_state


@fixture
def betse_sim_conf(
    betse_sim_conf_default: SimConfTestInternal) -> SimConfTestInternal:
    '''
    Per-test fixture creating a temporary minified simulation configuration
    file and returning a wrapper around this file.

    Configuration Modifications (On-disk)
    ----------
    This fixture copies the default simulation configuration file for this
    application, complete with all external assets (e.g., geometry masks)
    referenced and required by this file, into a temporary directory whose
    basename is the name of the test requesting this fixture excluding the
    prefixing substring ``test_``. When requested by the
    ``test_cli_sim_default`` test, for example, this fixture creates a
    temporary simulation configuration file
    ``{tmpdir}/cli_sim_default/sim_config.yaml`` for the absolute path
    ``{tmpdir}`` of this test session's root temporary directory (e.g.,
    ``/tmp/pytest-0/cli_sim_default/sim_config.yaml``).

    This directory and thus simulation configuration is safely accessible
    *only* for the duration of the current test. Subsequently run tests and
    fixtures *cannot* safely reuse this configuration.

    Configuration Modifications (In-memory)
    ----------
    This fixture also transforms the in-memory instance of the
    :class:`betse.science.parameters.Parameters` class encapsulating this
    configuration as follows:

    * All configuration options either requiring interactive input *or*
      displaying interactive output are disabled (e.g., plots, animations).
    * The space and time costs associated with simulating this configuration
      are safely minimized in a manner preserving all features.

    Since this fixture does *not* write these changes back to this file, the
    parent fixture or test is expected to do so manually (e.g., by calling the
    :meth:`SimConfTestInternal.p.save_inplace` method on the object returned
    by this fixture).

    Parameters
    ----------
    betse_sim_conf_default : SimConfTestInternal
        Wrapper around a temporary non-minified simulation configuration file.

    Returns
    ----------
    SimConfTestInternal
        Wrapper around a temporary simulation configuration
        file specific to the current test, including such metadata as:

        * The absolute path of this configuration's on-disk YAML file.
        * This configuration's in-memory dictionary deserialized from this
          file.
    '''

    # Minimize the space and time costs associated with this configuration.
    betse_sim_conf_default.config.minify()

    # Return this wrapper.
    return betse_sim_conf_default


@fixture
def betse_sim_conf_compat(
    betse_temp_dir: LocalPath) -> SimConfTestInternal:
    '''
    Per-test fixture creating and returning a wrapper around a temporary
    simulation configuration file (complete with a pickled seed,
    initialization, and simulation) produced by the oldest version of this
    application for which the current version of this application guarantees
    backward compatibility.

    Caveats
    ----------
    Unlike the object returned by the comparable :func:`betse_sim_conf`
    fixture, the object returned by this fixture is *not* safely modifiable by
    the current version of this application. Doing so would invalidate the
    pickled files produced by the older version of this application, which
    would largely defeat the purpose of invoking this fixture.

    Parameters
    ----------
    betse_temp_dir : LocalPath
        Wrapper around a temporary directory isolated to the current test.

    Returns
    ----------
    SimConfTestInternal
        Wrapper around a temporary simulation configuration file specific to
        the current test, complete with pickled seed, initialization, and
        simulation files produced by the older version of this application.
    '''

    # Defer heavyweight imports.
    from betse import metadata
    from betse.util.path import files
    from betse.util.py.module import pymodule

    # Initialize the application metadata singleton, as required by the
    # subsequent instantiation of the "SimConfTestInternal" subclass.
    initter.init_app()

    # Absolute dirname of the directory of the root "betse_test" package.
    betse_test_dirname = pymodule.get_dirname('betse_test')

    # Absolute filename of the default simulation configuration file produced
    # by the oldest version of this application for which the current version
    # of this application guarantees backward compatibility.
    src_conf_filename = files.join_or_die(
        betse_test_dirname,
        'data',
        metadata.GIT_TAG_COMPAT_OLDEST,
        'yaml',
        'sim_config.yaml',
    )

    # Wrapper wrapping the default simulation configuration file copied into
    # this temporary directory and sanitized therein.
    sim_state = SimConfTestInternal(
        src_conf_filename=src_conf_filename,
        trg_conf_filepath=betse_temp_dir.join('sim_config.yaml'))

    # Return this wrapper *WITHOUT* minifying this configuration.
    return sim_state
