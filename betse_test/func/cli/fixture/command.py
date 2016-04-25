#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
External command fixtures.

These fixtures automate testing of the current versions of BETSE's externally
runnable CLI and GUI commands (e.g., `betse`, `betse-qt4`) in subprocesses of
the current test process, regardless of whether these commands have been
editably installed (i.e., as synchronized symlinks rather than desynchronized
copies) into the current Python environment or not.
'''

# ....................{ IMPORTS                            }....................
from pytest import fixture

# ....................{ CLASSES                            }....................
class CLIRunner(object):
    '''
    CLI-specific test runner, efficiently executing the external command for the
    BETSE CLI (i.e., `betse`) in the active Python interpreter.

    For both efficiency and reliably, this runner does _not_ actually execute
    this command. Doing so would introduce installation complications and
    portability concerns (e.g., conflicting versions of the `betse` command in
    the current `${PATH}`). Instead, this runner:

    * Imports the `betse.cli.__main__` module implementing the BETSE CLI.
    * Passes this module's `run()` method the passed arguments.

    CLI-specific fixtures typically return instances of this class as a
    means of communicating this runner to other fixtures and tests.
    '''

    #FIXME: Implement me to accept zero or more arguments, import
    #"betse.cli.__main__", and pass these arguments to that entry point.

    def run(*args) -> None:
        '''
        Call the entry point for the BETSE CLI with the passed arguments.

        Parameters
        ----------
        args : list
            List of zero or more arguments to be passed to this entry point,
            corresponding exactly to the set of command-line arguments accepted
            by the external command for the BETSE CLI (i.e., `betse`).

        See Also
        ----------
        `betse --help`
            Further details on such arguments.
        '''

        pass

# ....................{ FIXTURES ~ low-level               }....................
@fixture(scope='session')
def betse_cli(tmpdir_factory, request) -> CLIRunner:
    '''
    Fixture returning a singleton instance of the `CLIRunner` class.

    Returns
    ----------
    CLIRunner
        This instance.
    '''

    return CLIRunner()
