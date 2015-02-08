#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

#FIXME; The following snippet courtesy Matthew Leingan affords an elegant means
#of integrating built-in Python argument parsing and logging;
#
#    import argparse
#    import logging
#
#    parser = argparse.ArgumentParser()
#    parser.add_argument('-d','--debug',
#        help='Print lots of debugging statements',
#        action="store_const",dest="loglevel",const=logging.DEBUG,
#        default=logging.WARNING
#    )
#    parser.add_argument('-v','--verbose',
#        help='Be verbose',
#        action="store_const",dest="loglevel",const=logging.INFO
#    )
#    args = parser.parse_args()
#    logging.basicConfig(level=args.loglevel)
#
#In his own words: "So if --debug is set, the logging level is set to DEBUG. If
#--verbose, logging is set to INFO. If neither, the lack of --debug sets the
#logging level to the default of WARNING."
#
#That said, devoting a separate command-line option to "--debug" strikes us as
#overkill when we can simply repeat "-vv" to achive the same, ala this snippet;
#
#    # This isn't quite right, as it uses "--verbose=1"-style integer
#    # assignment. We just want to repeat "-v" and/or "--verbose". Much simpler.
#    # Nonetheless, this should be of some use.
#    parser = argparse.ArgumentParser()
#    parser.add_argument("-v", "--verbose", const=1, default=0, type=int, nargs="?",
#                        help="increase verbosity: 0 = only warnings, 1 = info, 2 = debug. No number means info. Default is no verbosity.")
#    args = parser.parse_args()
#
#    logger = logging.getLogger()
#    if args.verbose == 0:
#        logger.setLevel(logging.WARN)
#    elif args.verbose == 1:
#        logger.setLevel(logging.INFO)
#    elif args.verbose == 2:
#        logger.setLevel(logging.DEBUG)

'''`betse`'s command line interface (CLI).'''

# ....................{ IMPORTS                            }....................
from betse.cli.cli import CLI
import argparse

# ....................{ MAIN                               }....................
def main() -> int:
    '''Run `betse`'s command line interface (CLI).

    This function is provided as a convenience to callers requiring procedural
    functions rather than classical methods (e.g., `setuptools`).

    Returns
    ----------
    int
        Exit status of such interface. This is a non-negative integer in
        `[0, 255]` where 0 signifies success and all other values failure.
    '''
    return CLICLI().run()

# ....................{ MAIN                               }....................
class CLICLI(CLI):
    '''`betse`'s command line interface (CLI).

    Attributes
    ----------
    '''
    def __init__(self):
        pass

    def _run(self) -> None:
        '''
        Run `betse`'s command line interface (CLI).
        '''
        # Parse command-line arguments into object attributes.
        self._parse_args()
        # self._parse_common_args()

        #FIXME: Implement me.
        # Run the command specified by such arguments.

    def _parse_args(self):
        '''Parse all currently passed command-line arguments.
        '''
        #FIXME: Display a default help message, when the user passes no
        #arguments.

        # Define global command-line arguments (i.e., arguments applicable to
        # all subparsers, added below).
        parser = argparse.ArgumentParser(
            # Program name and description.
            prog = ''.join((metadata.NAME, ' ', metadata.__version__,)),
            description = metadata.DESCRIPTION,
            usage = ''.join((
                metadata.SCRIPT_NAME_CLI, ' <command> [<arg>...]',
            ))
        )
        self._add_parser_common_args(parser)

        # Define argument subparsers (i.e., actions).
        subparsers = parser.add_subparsers(
            title = 'actions',
            description = 'Actions to be performed.',
        )

        #FIXME: Not entirely clear as to why we require or want this dictionary.

        # Dictionary from action name to subparser object.
        action_parsers = {}
#        action_parsers['help'   ] = self._add_subparser_help(subparsers)
        action_parsers['run'] = self._add_subparser_run(subparsers)

        # Add an identifying name and description for each action parser.
        for action_parser_name, action_parser in action_parsers.iteritems():
            assert action_parser, '"{}" is None'.format(action_parser_name)
            action_parser.set_defaults(action_name=action_parser_name)

        #FIXME: This doesn't seem quite right. Don't we want to parse arguments
        #only after defining all argument syntax for actions as well?

        # Parse command-line arguments according to such specification.
        parser.parse_args(namespace=self)

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
