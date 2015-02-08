#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

#FIXME: Improve output on uncatched exceptions by wrapping all CLI operations in
#a try-except block that:
#
#* Catches all exceptions.
#* Improves such output.
#* Exits the current process with status 1.

'''Abstract command line interface (CLI).'''

# ....................{ IMPORTS                            }....................
from abc import ABCMeta
from betse import dependencies
from betse.util.path import dirs, log
from betse.util import output
import argparse, traceback

# ....................{ MAIN                               }....................
class CLI(metaclass = ABCMeta):
    '''
    Abstract command line interface (CLI) suitable for use by both the front-
    facing CLI and GUI applications for `betse`.

    Attributes
    ----------
    _logger : Logger
        Low-level object implementing such logging.
    '''
    def __init__(self):
        #FIXME: Initialize me, please.
        self._logger = None

    def run(self) -> int:
        '''
        Run the command line interface (CLI) defined by the current subclass.

        Returns
        ----------
        int
            Exit status of such interface, guaranteed to be a non-negative
            integer in `[0, 255]`, where 0 signifies success and all other
            values failure.
        '''
        try:
            # If at least one mandatory dependency is missing, fail.
            dependencies.die_unless_satisfied()

            # Make betse's top-level dot directory if not found.
            dirs.make_unless_found(dirs.DOT_DIR)

            # Parse command-line arguments and run the specified command.
            self._run()

            # Exit with successful exit status from the current process.
            return 0
        except Exception as exception:
            # Print such exception.
            self._print_exception(exception)

            # Exit with failure exit status from the current process. If such
            # exception provides a system-specific exit status, use such status;
            # else, use the default such status.
            return getattr(exception, 'errno', 1)

    def _print_exception(self, exception: Exception) -> None:
        '''
        Print the passed exception to standard error *and* log such exception.
        '''
        # If such printing itself raises an exception...
        try:
            assert isinstance(exception, Exception),\
                '"{}" not an exception.'.format(exception)

            #FIXME: Print such exception.

            # If a logger has been initialized, log such exception *AFTER*
            # printing such exception to standard error. (Logging is
            # substantially more fragile and hence likely to itself raise
            # further exceptions.)
            if self._logger:
                self._logger.exception(exception)
        # ...catch and print such exception using standard Python facilities
        # guaranteed not to raise additional exceptions.
        except Exception:
            output.error('print_exception() recursively raised exception:\n')
            traceback.print_exc()

    @abc.abstractmethod
    def _run(self):
        '''
        Parse all passed command-line arguments and run the specified command.
        '''
        pass

# --------------------( WASTELANDS                         )--------------------
        # return exception.errno if hasattr(exception, 'errno') else 1
#FUXME: Define a new module "betse/dependency.py" performing validation of
#external dependencies, both Python and non-Python. Although we believe "yppy"
#implemented such functionality, google about for the optimum Python 3 solution
#to this presumably commonplace problem.

    # _args : list
    #     List of zero or more arguments passed to such interface (e.g., from the
    #     command line).
        # Initialize such arguments to the current argument list, excluding
        # such list's first item. By cross-platform consent, such item is
        # *ALWAYS* the command name for the current process (e.g., "betse") and
        # hence ignorable.
        # self._args = sys.argv[1:]

#List of zero or more external arguments passed from the command line.
# def main(args = None):
#     CLI().run()
#     '''Run betse`'s command line interface (CLI).
#
#     Parameters
#     ----------
#     args : list, optional
#         List of zero or more arguments passed to such interface (e.g., from the
#         command line) or `None` if called as the entry point in an external
#         script installed by `setuptools`.
#     '''
#     # If called from a setuptools-installed script, copy such arguments from the
#     # argument list excluding the first item of such list. By cross-platform
#     # agreement, such item is *ALWAYS* the command name of the current process
#     # (e.g., "betse") and hence ignorable.
#     if args is None:
#         args = sys.argv[1:]

# if __name__ == '__main__':
#     main()

#FUXME; Configure me for CLI usage. Note that I'm no longer convinced that the
#way we launched "yppy" (e.g., "bin/yppy.bash") was ideal. We really want to do
#the "Pythonic" thing here. ruamel.yaml, for example, installs a Python wrapper
#"/usr/lib/yaml" which (in order):
#
#* Finds an appropriate Python interpreter.
#* Replaces the current process with the result of interpreting
#  "/usr/lib/python-exec/python${PYTHON_VERSION}/yaml". Such file appears to be
#  autogenerated by setuptools at installation time.
#FUXME; Hmm: it looks like we want a new file "betse/__main__.py" resembling:
#    from betse.main import main
#    main()
#This then permits betse to be run as follows:
#    # Yes, you either have to be in the parent directory of the directory
#    # containing such "__main__.py" file *OR* you have to fiddle with
#    # ${PYTHONPATH}.
#    >>> cd ~/py/betse
#    >>> python -m betse
#Naturally, this lends itself well to shell scripting. (Yay!)
#FUXME; Wo! Even nicer. setuptools has implicit support for "__main__.py"-style
#entry points. We just need a "setup.py" resembling:
#    setup(
#        # [...]
#        entry_points={
#            'betse': ['betse = betse.main:main'],
#        },
#    )
#What's sweet about this is that we can define additional separate scripts with
#deeper entry points if we need and or want to.
