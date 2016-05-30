#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Functional tests for BETSE's CLI testing all simulation-specific subcommands
(e.g., `betse try`).
'''

# ....................{ IMPORTS                            }....................
# from betse_test.mark.skip import skip

# ....................{ TESTS                              }....................
#FIXME: While useful, this should be split apart into its constituent steps,
#most of which should probably become fixtures for efficient reuse. Wait, no.
#Ideally, we should construct a linear chain of tests, each exercising a
#simulation phase (e.g., "seed", "init", "sim") depending on the success of the
#prior test in this chain.
#FIXME: Right. No. This is clearly a job for parametrization. Ideally, we only
#want to declare a single test for each simulation configuration which accepts
#the names of all CLI subcommands to be passed the configuration file. Perhaps
#we want to replace the "betse_cli" fixture requested below with a new
#"betse_cli_sim" fixture producing a class inheriting the class produced by the
#former and adding support for parametrization? Consider it.
#
#Indeed, that's the way. However, we don't want tests themselves to specify
#parameters. Instead, we want the new "betse_cli_sim" fixture to unconditionally
#parametrize *ITSELF* by the in-order list of all CLI subcommands to be run:
#e.g.,
#
#    (
#        ('seed'), ('init'), ('sim'),
#        ('plot', 'seed'), ('plot', 'init'), ('plot', 'sim'),
#    )
#
#Why? A couple obvious reasons:
#
#1. We don't want to duplicate that list for every test. By centralizing that
#   list into a single fixture, changes are dramatically simplified.
#2. Tests may *THINK* they know which CLI subcommands they need run, but...
#   really, they don't. Each new configuration really requires everything to be
#   revalidated from scratch. This approach forces that.
#3. By parametrizing the "betse_cli_sim" fixture, all tests requesting that
#   fixture will themselves receive the same parameters -- which means that
#   users can interactively test exact CLI subcommands for exact configurations
#   by simply using "./test -k". Great!
#
#In short, this is the way forward.
#FIXME: There are numerous issues with the parametrization approach, including:
#* Interpameter dependencies. How do we ensure that a prior parameter is run
#  before the current parameter to be run? For example, if the user runs "./test
#  -k test_cli_sim_default[plot_sim]", then how do we ensure that the "sim"
#  parameter is run first? While it probably would be feasible to extend hpk42's
#  class-based approach to parameters, it hardly seems worthwhile.
#* Phase xfailing. Under the parameter approach, how would be mark an individual
#  phase of a simulation configuration -- say, "plot_sim", as xfailing? We have
#  no idea.
#
#In short, a class approach appears to make considerably more sense.
#FIXME: Let's a-go:
#
#* Define a new "CLISimTesterABC" base class inheriting
#  "betse_test.util.testabc.SerialTestABC".
#* In "CLISimTesterABC", define one method for each prospective phase requiring
#  a "betse_cli" fixture: e.g.,
#      def test_cli_sim_plot_sim(self, betse_cli) -> None:
#  Each such method should access the new "self.sim_state" attribute. See below.
#* Set "CLISimTesterABC.__test__ = False" to force these methods to be ignored.
#* Define a new "CLISimDefaultTester" subclass inheriting "CLISimTesterABC".
#* Set "CLISimDefaultTester.__test__ = True" to unignore these methods.
#* The only question then becomes: what becomes of the
#  "betse_sim_config_default" fixture? The simple solution might be to:
#  * Refactor "betse_sim_config_default" into an autouse-style class fixture
#    defined in "CLISimDefaultTester". On execution, this fixture should add a
#    new "cls.sim_state" attribute refering to the returned fixture object. In
#    theory, test methods may then access this attribute.
#  * The only question then becomes: how does "CLITestRunner" inspect for the
#    "betse_sim_config_default" fixture? Specifically, are fixtures that are
#    implicitly run as autouse added to the "request.fixturenames" list? We
#    strongly suspect the answer is *YES*, thank Odin.
#FIXME: This is a fairly heavyweight approach. From that perspective, the
#parameters- based approach would certainly be preferable. My reasoning for
#preferring the class- to parameters-based approach were fairly flimsy, frankly.
#The only question then becomes: can we refactor the existing
#"betse_test.conftest" hooks to support parameters-based serial testing?
#FIXME: Automatically coerce the "world size" configuration option to "75e-6".

def test_cli_sim_default(
    betse_cli, betse_sim_config_default) -> None:
    # betse_cli_sim, betse_sim_config_default) -> None:
    '''
    Test the simulation of the default simulation configuration.

    Parameters
    ----------
    betse_cli_sim : CLITestRunner
        Test-specific object encapsulating the simulation-specific BETSE CLI.
    betse_sim_config : SimTestState
        Test-specific object encapsulating this simulation configuration file.
    '''

    # Pass the basename rather than absolute path of this configuration file,
    # testing that the "betse_sim_config_default" fixture has set the current
    # working directory (CWD) to the directory containing this file.
    betse_cli('sim', betse_sim_config_default.config.basename)
    # betse_cli_sim.run()
    pass
