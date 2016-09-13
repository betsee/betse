#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Functional tests for BETSE's CLI testing all simulation-specific subcommands
(e.g., `betse try`).
'''

# ....................{ IMPORTS                            }....................
from betse_test.util.mark.param import parametrize_test

# ....................{ TESTS                              }....................
#FIXME: Refactor the betse_cli_sim() fixture to be:
#
#* Indirectly rather than directly parametrized by the currently accepted
#  tuple of subcommand shell words.
#FIXME: To do so sanely, however, we'll *ALSO* need to dynamically synthesize
#the passed list of parametrization IDs based on the passed shell words. That,
#in turn, would appear to necessitate the definition of a new
#@parametrize_test_cli_sim_config decorator based on the existing
#@parametrize_test decorator. The former should wrap the latter with additional
#functionality for:
#
#* Accepting the following high-level keyword arguments:
#  * "subcommand_shell_words", the mandatory tuple of tuples of subcommand
#    shell words to run.
#  * "sim_config_modifier", the optional callable modifying the passed
#    simulation configuration object.
#* Synthesizing the following lower-level keyword arguments:
#  * "ids" from the contents of "subcommand_shell_words".
#  * "fixtures" to parametrize the requisite:
#    * "betse_sim_config" fixture from the contents of "sim_config_modifier".
#      Note that this fixture *MUST* be parametrized the same number of times as
#      the "betse_cli_sim" fixture -- which probably means repeating the
#      "sim_config_modifier" callable len(subcommand_shell_words) times as
#      elements of a tuple. Trivial... but let's not forget!
#    * "betse_cli_sim" fixture from the contents of
#      "subcommand_shell_words".
#FIXME: For convenience, we might also consider forcefully prepending the
#keyword arguments accepted by this test by "betse_cli_sim" and
#"betse_sim_config". It's unclear whether or not doing so would actually be
#acknowledged by whatever introspection py.test uses to discover fixtures
#requested by tests, however. In short, this is hardly essential and probably
#won't work as desired.
#FIXME: O.K.; sadly, all of the above is off the table. Py.test appears to be
#fundamentally incapable of combining multiple fixtures of different scope
#accepting different indirect parametrizations. While there probably is a means
#of eventually doing so, several observations are beginning to become apparent:
#
#* The current disconnection between the "betse_cli_sim" fixture and
#  ""betse_sim_config_*" family of fixtures is, frankly, absurd. The former goes
#  to excrutiating lengths to access the latter.
#* The difficulty of indirectly parametrizing fixtures, even given our new
#  @parametrize_test() fixture.
#* The current infeasibility of reliably serializing parametrized fixtures and
#  the inability to use the "xdist" plugin that results from that infeasibility.
#
#All of this suggests, sadly, that our design fundamentally went down the wrong
#road and no longer scales at all. Our fault. Basically, we need to
#*DRAMATICALLY* simplify this entire design as follows:
#
#* Cease attempting to run BETSE subcommand phases as parametrized serial
#  independent tests.
#* Instead, define a new "CLISimTester" class loosely based on the existing
#  "CLITesterPreArged" class as follows:
#  * Generalize the CLISimTester.__init__() method to accept no
#    "subcommand_args" parameter.
#  * Shift the contents of the SimTestState.__init__() and
#    SimTestState.get_command_context() methods to the beginning of the
#    CLISimTester.run() method. To do so safely, we'll want to embed the call to
#    dirs.current() in a "with"-style context manager. Trivial!
#  * Rename the CLISimTester.run() method to.run_subcommands().
#  * Refactor this method to accept the following keyword arguments:
#    * "subcommands_args", a mandatory sequence of sequences of shell words to
#      run.
#    * "config_modifier", an optional lambda modifying this simulation
#      configuration.
#    For each such sequence of shell words, this method should:
#    * Log an "info" message declaring the subcommand to be run.
#    * Run that subcommand.
#  * Define a new CLISimTester.run_subcommands_default() method to accept only
#    the optional "config_modifier" argument, internally deferring to the
#    run_subcommands() method by defaulting the "subcommands_args" argument to
#    the tuple of all default subcommands.
#* Refactor the "betse_cli_sim" fixture to return instances of this new class.
#* Refactor all tests below to:
#  * Request *ONLY* the "betse_cli_sim" fixture.
#  * Call betse_cli_sim.run_subcommands(), optionally passed a "config_modifier"
#    lambda.
#* Dramatically simplify the "betse_cli" fixture and "CLITester" class by
#  removing all mention of "contexts". This was a terrible idea. That should
#  have been immediately apparent. Now it is.
#* Remove all "betse_sim_config*" fixtures.
#* Remove the now-obsoleted "SimTestState" and "CLITesterPreArged" classes.
#* Re-enable "xdist" support. Maybe? If we do so, we should only do so under the
#  caveat that "xdist" support is disabled when "-s" is passed. (Pretty sure we
#  currently ensure this, but... let's make sure.)

def test_cli_sim_default(
    betse_cli_sim: 'CLITesterPreArged',
    betse_sim_config_default: 'SimTestState',
) -> None:
    '''
    Test simulating the default simulation configuration.

    Parameters
    ----------
    betse_cli_sim : CLITesterPreArged
        Object encapsulating CLI-driven simulation testing.
    betse_sim_config_default : SimTestState
        Object encapsulating the default simulation configuration.
    '''

    # Test the currently parametrized simulation-specific BETSE CLI subcommand.
    betse_cli_sim()


def test_cli_sim_visuals(
    betse_cli_sim: 'CLITesterPreArged',
    betse_sim_config_visuals: 'SimTestState',
) -> None:
    '''
    Test simulating all exported visuals (e.g., plots, animations) _and_
    simulation features required by these visuals.

    Parameters
    ----------
    betse_cli_sim : CLITesterPreArged
        Object encapsulating CLI-driven simulation testing.
    betse_sim_config_visuals : SimTestState
        Object encapsulating the simulation configuration enabling these plots
        and animations.
    '''

    # Test the currently parametrized simulation-specific BETSE CLI subcommand.
    betse_cli_sim()
