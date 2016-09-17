#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Functional tests for BETSE's CLI testing all simulation-specific subcommands
(e.g., `betse try`).
'''

# ....................{ IMPORTS                            }....................
import pytest
from betse_test.util.mark.skip import skip_unless_matplotlib_anim_writer
from betse_test.util.mark.fail import xfail

# ....................{ TESTS                              }....................
def test_cli_sim_default(betse_cli_sim: 'CLISimTester') -> None:
    '''
    Test simulating the default simulation configuration.

    Parameters
    ----------
    betse_cli_sim : CLISimTester
        Object running multiple BETSE CLI simulation subcommands.
    '''

    # Persist all in-memory configuration changes back to disk.
    betse_cli_sim.sim_state.config.overwrite()

    # Test all default simulation-specific subcommands with this configuration.
    betse_cli_sim.run_subcommands_default()


def test_cli_sim_visuals(betse_cli_sim: 'CLISimTester') -> None:
    '''
    Test simulating all exported visuals (e.g., plots, animations) _and_
    simulation features required by these visuals.

    Parameters
    ----------
    betse_cli_sim : CLISimTester
        Object running multiple BETSE CLI simulation subcommands.
    '''

    # Enable all exported visuals and features required by these visuals.
    betse_cli_sim.sim_state.config.enable_visuals_all()

    # Persist all in-memory configuration changes back to disk.
    betse_cli_sim.sim_state.config.overwrite()
    # print('\n!!!!!after solving: {}'.format(betse_cli_sim.sim_state.config._config['results options']['after solving']))

    # Test all default simulation-specific subcommands with this configuration.
    betse_cli_sim.run_subcommands_default()


# Sadly, all existing higher-level parametrization decorators defined by the
# "betse_test.util.mark.params" submodule fail to support embedded py.test-style
# "skipif" and "xfail" markers. Consequently, we leverage the lower-level
# parametrization decorator shipped with py.test itself.
@pytest.mark.parametrize(
    ('writer_name', 'filetype'), (
        skip_unless_matplotlib_anim_writer('avconv')(('avconv', 'mp4')),
        skip_unless_matplotlib_anim_writer('ffmpeg')(('ffmpeg', 'mkv')),

        #FIXME: Research this deeper, please. Are all Mencoder-based writers
        #genuinely broken (doubtful), is this our fault (very possible), or is
        #this a platform- or version-specific issue and hence mostly not our
        #fault (also very possible)?
        # skip_unless_matplotlib_anim_writer('mencoder')(('mencoder', 'avi')),
        xfail(reason='Mencoder-based writers fail with obscure errors.')(
            ('mencoder', 'avi')),

        # ImageMagick only supports encoding animated GIFs and hence is
        # effectively a joke writer. Since it remains supported, however, we
        # test it with a reasonable facsimile of a straight face.
        skip_unless_matplotlib_anim_writer('imagemagick')(
            ('imagemagick', 'gif')),
    ),
)
def test_cli_sim_video(
    betse_cli_sim: 'CLISimTester', writer_name: str, filetype: str) -> None:
    '''
    Test simulating at least one animation (and all simulation features required
    by that animation) and encoding that animation as compressed video of the
    parametrized filetype via the matplotlib animation writer of the
    parametrized name.

    Since video encoding is justifiably expensive in both space and time, this
    test encodes only a single animation (rather than all available animations)
    for a single simulation-specific BETSE CLI subcommand (rather than all
    available subcommands). The existing :func:`test_cli_sim_visuals` test
    already exercises all animations for all subcommands, albeit with video
    encoding disabled.

    Parameters
    ----------
    betse_cli_sim : CLISimTester
        Object running multiple BETSE CLI simulation subcommands.
    writer_name : str
        Name of the matplotlib animation writer with which to encode video
        (e.g., `ffmpeg`, `imagemagick`).
    filetype : str
        Filetype of videos to encode with this writer (e.g., `mkv`, `mp4`).
    '''

    # Enable encoding of one animation to this filetype with this writer..
    betse_cli_sim.sim_state.config.enable_anim_video(writer_name, filetype)

    # Persist all in-memory configuration changes back to disk.
    betse_cli_sim.sim_state.config.overwrite()

    # Test the minimum number of simulation-specific subcommands required to
    # exercise video encoding with this configuration. Since the "init"
    # subcommand requiring the "seed" subcommand satisfies this constraint, all
    # subsequent subcommands (e.g., "sim", "plot init") are omitted.
    betse_cli_sim.run_subcommands(('seed',), ('init',),)
