#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
CLI-specific functional tests exercising **solver-agnostic simulations** (i.e.,
simulations arbitrarily supporting all simulation solvers, including both the
complete BETSE solver *and* the "fast" equivalent circuit solver).
'''

# ....................{ IMPORTS                           }....................
import pytest
from betse.util.test.pytest.mark.pytfail import xfail
from betse.util.test.pytest.mark.pytskip import (
    skip_unless_matplotlib_anim_writer, skip_if_requirement)

# ....................{ TESTS                             }....................
#FIXME: Sadly, our current approach to backward compatibility testing is
#fundamentally flawed. Why? Because third-party dependencies (e.g., pytest,
#setuptools) continue to break backward compatibility. Ironically, this renders
#our own attempts to preserve backward compatibility infeasible.
#
#Technically, we *COULD* probably circumvent this issue by installing the older
#version of BETSE checked out for this functional test within a virtual
#environment of some sort (e.g., Pipenv). Doing so would invest even more
#scarce development resources in a probably flawed testing regime, however.
#
#Pragmatically, the optimal approach is simply to embed a copy of the
#"betse.data" subpackage corresponding to that of the oldest version of BETSE
#with which we preserve backward compatibility within a new "betse_test.data"
#subpackage -- presumably gated by version-specific subdirectories: e.g.,
#
#* "betse_test.data.0_5_2", containing the exact subset of files provided by
#  the "betse.data.yaml" subdirectory of BETSE 0.5.2 required to reproduce
#  this test's requirements.
#* "betse_test.data.0_6_0", likewise for BETSE 0.6.0.
#
#In other words: *WHAT WERE WE THINKING.* Well, O.K.; we knew what we were
#thinking. We were attempting to avoid data duplication. In this case, the
#minimal set of all data required to safeguard backward compatibility is
#probably of ignorable filesize. In simpler words, we chose poorly.

# This function test is well-known to be incompatible with recent versions of
# py.test, raising exceptions resembling:
#
#     ============================= test session starts ==============================
#     platform linux -- Python 3.7.1, pytest-4.1.1, py-1.7.0, pluggy-0.8.0 -- /builds/betse/betse/conda-env/bin/python
#     cachedir: .pytest_cache
#     rootdir: /tmp/pytest-of-root/pytest-0/cli_sim_compat0/betse_old, inifile: pytest.ini
#     collecting ...
#     ==================================== ERRORS ====================================
#     _________________ ERROR collecting betse_test/func/test_cli.py _________________
#     /builds/betse/betse/conda-env/lib/python3.7/site-packages/pluggy/hooks.py:284: in __call__
#         return self._hookexec(self, self.get_hookimpls(), kwargs)
#     /builds/betse/betse/conda-env/lib/python3.7/site-packages/pluggy/manager.py:67: in _hookexec
#         return self._inner_hookexec(hook, methods, kwargs)
#     /builds/betse/betse/conda-env/lib/python3.7/site-packages/pluggy/manager.py:61: in <lambda>
#         firstresult=hook.spec.opts.get("firstresult") if hook.spec else False,
#     /builds/betse/betse/conda-env/lib/python3.7/site-packages/_pytest/python.py:225: in pytest_pycollect_makeitem
#         res = list(collector._genfunctions(name, obj))
#     /builds/betse/betse/conda-env/lib/python3.7/site-packages/_pytest/python.py:405: in _genfunctions
#         self.ihook.pytest_generate_tests(metafunc=metafunc)
#     /builds/betse/betse/conda-env/lib/python3.7/site-packages/pluggy/hooks.py:284: in __call__
#         return self._hookexec(self, self.get_hookimpls(), kwargs)
#     /builds/betse/betse/conda-env/lib/python3.7/site-packages/pluggy/manager.py:67: in _hookexec
#         return self._inner_hookexec(hook, methods, kwargs)
#     /builds/betse/betse/conda-env/lib/python3.7/site-packages/pluggy/manager.py:61: in <lambda>
#         firstresult=hook.spec.opts.get("firstresult") if hook.spec else False,
#     /builds/betse/betse/conda-env/lib/python3.7/site-packages/_pytest/python.py:132: in pytest_generate_tests
#         metafunc.parametrize(*marker.args, **marker.kwargs)
#     /builds/betse/betse/conda-env/lib/python3.7/site-packages/_pytest/python.py:892: in parametrize
#         function_definition=self.definition,
#     /builds/betse/betse/conda-env/lib/python3.7/site-packages/_pytest/mark/structures.py:114: in _for_parametrize
#         if len(param.values) != len(argnames):
#     E   TypeError: object of type 'MarkDecorator' has no len()
#
# Since ours appears to be the only py.test-based test suite exhibiting this
# exception, identifying and resolving the underlying culprit (e.g., by
# monkey-patching) is effectively infeasible. Moreover, since it remains
# unclear which py.test version introduced this incompatibility, we have little
# choice but to skip the entire py.test 4.x release line and hope for the best.
#FIXME: After simplifying this functional test as detailed above, remove all of
#the following:
#
#* The "--export-sim-conf-dir" option, defined by the pytest_addoption() hook
#  in the top-level "hetse_test.conftest" plugin.
#* The "test_sim_export" submodule.
@skip_if_requirement('pytest >= 4.0.0')
def test_cli_sim_compat( betse_cli_sim_compat: 'CLISimTester') -> None:
    '''
    Functional test exercising all simulation subcommands required to validate
    backward compatibility with a temporary simulation configuration file
    (complete with a pickled seed, initialization, and simulation) produced by
    the oldest version of this application for which the current version of
    this application guarantees backward compatibility.

    Design
    ----------
    Validating backward compatibility requires validating that the current
    version of this application can successfully load *all*:

    * Simulation configuration files loadable by this older version.
    * Pickled seeds, initializations, and simulations saved by this older
      version.

    Technically, running the:

    * ``plot sim`` subcommand would test the loadability of pickled
      simulations.
    * ``sim`` subcommand would test the loadability of pickled initializations.
    * ``init`` subcommand would test the loadability of pickled seeds.

    That said, running the ``init`` and ``sim`` subcommands would erroneously
    overwrite the pickled seeds and initializations saved by this older
    version. Doing so would invite edge-case issues that are best avoided; in
    particular, care would need to be taken to run the ``sim`` subcommand
    *after* the ``init`` subcommand. To ameliorate these concerns, plotting
    subcommands guaranteed to both test the loadability of pickled objects
    *and* be side effect-free are run instead.

    Parameters
    ----------
    betse_cli_sim_compat : CLISimTester
        Object running BETSE CLI simulation subcommands against a temporary
        simulation configuration produced by this older application version.
    '''

    # Test all simulation-specific plotting subcommands on this configuration.
    betse_cli_sim_compat.run_subcommands(
        *betse_cli_sim_compat.SUBCOMMANDS_PLOT,

        # Avoid overwriting the previously exported simulation configuration.
        is_overwrite_conf=False)

    #FIXME: We unavoidably break backward compatibility with respect to pickled
    #objects and hence currently only test backward compatibility with respect
    #to the simulation configuration. This is non-ideal, but sufficient for the
    #moment. Ideally, this should be reverted on releasing a new version.
    # betse_cli_sim_compat.run_subcommands(
    #     ('seed',),
    #     is_overwriting_config=False,)


def test_cli_sim_default(betse_cli_sim_default: 'CLISimTester') -> None:
    '''
    Functional test exercising the default simulation configuration on all
    simulation phases (i.e., seed, initialization, and simulation), principally
    for computational instability.

    Unlike the minified simulation configuration leveraged by all other tests,
    the default simulation configuration leveraged by this test is unmodified
    (except for disabling interactive features, which non-interactive testing
    unavoidably requires). In particular, this configuration is *not* minified
    for efficient usage by tests and thus incurs a significant performance
    penalty. This test exercises only the minimal subset of this configuration
    required to detect computational instabilities in this configuration --
    namely, the first two simulation phases.

    Parameters
    ----------
    betse_cli_sim_default : CLISimTester
        Object running non-minified BETSE CLI simulation subcommands.
    '''

    # Test only the first two simulation-specific subcommands with this
    # configuration, as documented above.
    betse_cli_sim_default.run_subcommands_sim()


# Sadly, all existing higher-level parametrization decorators defined by the
# "betse.util.test.pytest.mark.params" submodule fail to support embedded py.test
# "skipif" and "xfail" markers. Consequently, we leverage the lower-level
# parametrization decorator shipped with py.test itself.
@pytest.mark.parametrize(
    ('writer_name', 'filetype'), (
        pytest.param(
            'avconv', 'mp4',
            marks=skip_unless_matplotlib_anim_writer('avconv')),
        pytest.param(
            'ffmpeg', 'mkv',
            marks=skip_unless_matplotlib_anim_writer('ffmpeg')),

        #FIXME: Research this deeper, please. Are all Mencoder-based writers
        #genuinely broken (doubtful), is this our fault (very possible), or is
        #this a platform- or version-specific issue and hence mostly not our
        #fault (also very possible)?
        # skip_unless_matplotlib_anim_writer('mencoder')(('mencoder', 'avi')),
        pytest.param(
            'mencoder', 'avi', marks=xfail(
                reason='Mencoder-based writers fail with obscure errors.')),

        # ImageMagick only supports encoding animated GIFs and hence is
        # effectively a joke writer. Since it remains supported, however, we
        # test it with a reasonable facsimile of a straight face.
        pytest.param(
            'imagemagick', 'gif',
            marks=skip_unless_matplotlib_anim_writer('imagemagick')),
    ),
)
def test_cli_sim_video(
    betse_cli_sim: 'CLISimTester', writer_name: str, filetype: str) -> None:
    '''
    Functional test simulating at least one animation (and all simulation
    features required by that animation) and encoding that animation as
    compressed video of the parametrized filetype via the matplotlib animation
    writer of the parametrized name.

    Since video encoding is justifiably expensive in both space and time, this
    test encodes only a single animation (rather than all available animations)
    for a single simulation-specific BETSE CLI subcommand (rather than all
    available subcommands). The existing :func:`test_cli_sim_visuals` test
    already exercises all animations for all subcommands, albeit with video
    encoding disabled.

    Parameters
    ----------
    betse_cli_sim : CLISimTester
        Object running BETSE CLI simulation subcommands.
    writer_name : str
        Name of the matplotlib animation writer with which to encode video
        (e.g., `ffmpeg`, `imagemagick`).
    filetype : str
        Filetype of videos to encode with this writer (e.g., `mkv`, `mp4`).
    '''

    # Test the minimum number of simulation-specific subcommands required to
    # exercise video encoding with this configuration. Since the "init"
    # subcommand requiring the "seed" subcommand satisfies this constraint, all
    # subsequent subcommands (e.g., "sim", "plot init") are omitted.
    betse_cli_sim.run_subcommands(('seed',), ('init',),)
