#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2022 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Integration tests exercising the public API of the public
:mod:`betse.science.wrapper` subpackage.
'''

# ....................{ IMPORTS                            }....................
from betse_test._fixture.simconf.simconfclser import SimConfTestInternal
from py._path.local import LocalPath

# ....................{ TESTS                              }....................
def test_wrapper_default(betse_sim_conf: SimConfTestInternal) -> None:
    '''
    Integration test exercising default parameters accepted by the
    :meth:`betse.science.wrapper.BetseWrapper.__init__` method.

    Parameters
    ----------
    betse_sim_conf : SimConfTestInternal
        Object encapsulating a temporary simulation configuration file.
    '''

    # Defer test-specific imports.
    from betse.science.wrapper import BetseWrapper

    # # Save these changes back to the same file.
    # p.save_inplace()

    # BETSE wrapper configured by this file *AND* no other passed parameters,
    # intentionally exercising defaults set for those parameters.
    wrapper = BetseWrapper(config_filename=betse_sim_conf.p.conf_filename)

    # Run the default pipeline performing *ONLY* the initialization phase.
    wrapper.run_pipeline()


def test_wrapper_logging(
    betse_sim_conf: SimConfTestInternal,
    betse_temp_dir: LocalPath,
) -> None:
    '''
    Integration test exercising logging parameters accepted by the
    :meth:`betse.science.wrapper.BetseWrapper.__init__` method.

    Parameters
    ----------
    betse_sim_conf : SimConfTestInternal
        Object encapsulating a temporary simulation configuration file.
    betse_temp_dir : LocalPath
        Object encapsulating a temporary directory isolated to this test.
    '''

    # Defer test-specific imports.
    from betse.science.wrapper import BetseWrapper

    # Absolute filename of a logfile with arbitrary basename in this temporary
    # directory.
    log_file = betse_temp_dir.join('betse.log')

    # BETSE wrapper... *AND* no other passed parameters,
    # intentionally exercising defaults set for those parameters.
    wrapper = BetseWrapper(
        # Configured by this file.
        config_filename=betse_sim_conf.p.conf_filename,
        # Logging to another file.
        log_filename=str(log_file),
        # Logging *ALL* messages.
        log_level='ALL',
    )

    # Run the default pipeline performing *ONLY* the initialization phase.
    wrapper.run_pipeline()

    # Assert this logfile to have been created.
    assert log_file.check(file=1)

    # Assert this logfile to be non-empty.
    assert log_file.size() > 100
