#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Fixtures creating **shallow Git clones** (i.e., copies of a source Git
repository whose histories are truncated to the most recent commit to that
repository).
'''

#FIXME: Actually implement and document the fixture defined below. To do so,
#note that the following "git" command appears to implement our desired
#functionality:
#
#    git clone --branch v0.5.0 --depth 1 file:///home/leycec/py/betse betse_old
#
#Naturally, we'll want to:
#
#* Temporarily change the CWD to "str(betse_temp_dir)" *BEFORE* performing this
#  clone, as "git clone" only clones into the current directory.
#* Substitute the hard-coded source repository path given above with
#  pathtree.get_worktree_dirname_or_none(), which will require raising an
#  exception when that function returns "None".
#
#Note that the hard-coded target repository basename given above should be fine,
#as should literally anything, thanks to per-test fixture isolation.
#
#For generality, it might be wise to encapsulate this functionality into a new
#clone_git_worktree_shallow() function in a new "betse.util.path.gits"
#submodule. Overkill? Probably not. That function can perform useful checking to
#guarantee sanity. In particular, that function should:
#
#* Accept a target directory to perform this clone into and:
#  * Create this directory if needed.
#  * Temporarily change the CWD to this directory.
#
#It would probably then be useful to design an analogue to the
#simconfig/simconfer.betse_sim_config() fixture operating upon this... Hmmmmmmm.
#Complications arise, of course. Ideally, we would confine BETSE subcommands
#(e.g., "betse seed") to this clone for the duration of this test in a manner
#leveraging our existing test classes (e.g., "SimConfigTestWrapper"). Sadly,
#doing so sanely would require isolating this test to a "py.test" subprocess
#(which basically requires using "pytest-xdist", which is disfunctional for us)
#and munging "sys.path" to point to the clone.
#
#Since that won't work at all, we'll need to manually replicate the minimal
#amount of functionality performed by our test classes to get this to work
#*WITHOUT* reimplementing those classes in an insane manner. The only
#"SimConfigTestWrapper" method that appears essential for this use case is
#minify(), but we can probably do without even that (...at least initially).
#Err... no, we probably *CANNO* do without that. So, we'll want to do this
#manually by loading, minifying, and dumping the YAML file with our
#"betse.lib.yaml.yamls" helpers. Non-ideal, but such is life. *shrug*
#
#Let's just try doing this manually first by forking Python 3 subprocesses ala:
#
#    import os
#    from betse.util.path.command import cmdrun
#    from betse.util.py import pys
#
#    @fixture
#    def betse_cli_sim_git(...) -> None:
#
#        #FIXME: Close us up here by calling the
#        #clone_git_worktree_shallow() function delineated above.
#
#        #FIXME: Create a new simulation configuration from the default one by
#        #calling "betse config muh_sim.yaml" via the approach given below. Or,
#        #alternately, since this clone is per-test, we *COULD* probably just
#        #operate directly on "betse.data.yaml.sim_config.yaml". We would, of
#        #course, need to be *VERY* certain that we were operating on the
#        #correct "betse", which seems both fragile and dangerous. So, yeah:
#        #let's just call "betse config muh_sim.yaml" for safety. (Also,
#        #"betse.data.yaml.sim_config.yaml" is *NOT* intended to be used
#        #directly and may not necessarily behave as expected. So, don't.)
#
#        py_command_line_prefix = pys.get_command_line_prefix()
#
#        #FIXME: Parametrize this, of course! Wait! Maybe not. We *REALLY* don't
#        #want to perform this clone across multiple test parametrizations.
#        #Which implies we should just loop over "('seed', 'init', 'sim',)".
#        betse_subcommand = 'seed'
#        betse_command = py_command_line_prefix + (
#            '-m', 'betse.cli', betse_subcommand)
#
#        cmdrun.run_or_die(
#           command_words=betse_command,
#           popen_kwargs={
#               #FIXME: Specify an absolute rather than relative dirname here.
#               'env': dict(os.environ, PYTHONPATH='betse_old'),
#           },
#        )
#
#That looks fundamentally sane to me, but only ugly time will tell.

# ....................{ IMPORTS                            }....................
# from betse_test.util import requests
from pytest import fixture

# ....................{ FIXTURES                           }....................
# Test-scope fixture creating and returning a new object for each discrete test.
@fixture
def betse_git_clone(
    # request: '_pytest.python.FixtureRequest',
    betse_temp_dir: 'LocalPath',
) -> 'LocalPath':
    '''
    Per-test fixture creating a temporary directory isolated to the test
    requesting this fixture and returning an object encapsulating this
    directory.

    This directory is guaranteed to be specific to this test. The basename of
    this directory is the name of this test excluding the prefixing substring
    ``test_``. When requested by the ``test_cli_sim_default`` test, for example,
    this fixture creates a temporary directory ``{tmpdir}/cli_sim_default`` for
    the absolute path ``{tmpdir}`` of this test session's root temporary
    directory (e.g., ``/tmp/pytest-0/cli_sim_default``).

    This directory is safely accessible *only* for the duration of this test.
    Subsequently run tests and fixtures *cannot* safely reuse this directory,
    although doing so is technically feasible in unreliable edge-cases.

    Parameters
    ----------
    # request : _pytest.python.FixtureRequest
    #     Builtin fixture describing the parent fixture or test of this fixture.
    betse_temp_dir : LocalPath
        Object encapsulating a temporary directory isolated to the current test.

    Returns
    ----------
    LocalPath
        Object encapsulating this temporary directory.
    '''

    # Defer heavyweight imports.
    from betse.util.type.text import strs

    # Name of the current test.
    test_name = requests.get_tested_name(request)
    # print('    request.node: {}'.format(request.node))
    # print('    test_name: {}'.format(test_name))

    # Basename of the temporary directory containing this configuration file,
    # set to the name of the current test excluding the prefixing "test_".
    temp_dir_basename = strs.remove_prefix(
        text=test_name,
        prefix='test_',
        exception_message=(
            'Test name "{}" not prefixed by "test_".'.format(test_name)),
    )

    # Create this temporary directory and wrap this directory's absolute path
    # with a high-level "py.path.local" object. See also:
    #     http://pytest.org/latest/tmpdir.html#the-tmpdir-factory-fixture
    temp_dirpath = tmpdir_factory.mktemp(temp_dir_basename)

    # Return this object.
    return temp_dirpath
