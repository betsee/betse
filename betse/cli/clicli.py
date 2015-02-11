#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''`betse`'s command line interface (CLI).'''

# ....................{ IMPORTS                            }....................
from betse import metadata
from betse.cli.cli import CLI
from betse.util.io import stdout
from betse.util.path import files
from betse.util.system import processes
from collections import OrderedDict

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
        super().__init__()

    # ..................{ SUPERCLASS                         }..................
    def _configure_arg_parsing(self):
        #FIXME: Contemplate localizing.

        # Collection of argument subparsers parsing arguments for subcommands.
        self._arg_subparsers = self._arg_parser.add_subparsers(
            # Title of the subcommand section in help output.
            title = 'subcommands',

            # Description of the subcommand section in help output.
            description = 'Subcommand to be performed.',

            # Name of the attribute storing the passed command name.
            dest = 'command_name',
        )

        # Add argument subparsers accepting no options.
        self._arg_subparsers.add_subparser(
            name = 'info',
            description = 'print program metadata in key-value pair format',
        )

    def _run(self) -> None:
        '''
        Run `betse`'s command line interface (CLI).
        '''
        # If no subcommand was passed, print help output and return.
        if not self._args.command_name:
            self._arg_parser.print_help()
            return

        # Else, a subcommand was passed.
        #
        # Name of the method running such subcommand.
        subcommand_method_name = '_run_' + self._args.command_name

        # Method running such subcommand. If such method does *NOT* exist,
        # getattr() will raise a non-layman-readable exception. Typically, this
        # would be bad. In this case, however, argument parsing coupled with a
        # reliable class implementation guarantees such method to exist.
        subcommand_method = getattr(self, subcommand_method_name)

        # Run such subcommand.
        subcommand_method()

    # ..................{ SUBCOMMAND                         }..................
    def _run_info(self) -> None:
        '''
        Run the `info` subcommand.
        '''
        #FIXME; For aesthetics, convert to yppy-style "cli.memory_table" output.

        # Dictionary of string keys and string values to be output below,
        # ordered so as to preserve the specified order.
        info_key_to_value = OrderedDict((
            ('script name', processes.get_current_basename()),
            ('version',     metadata.__version__),
            ('config file', files.DEFAULT_CONFIG_FILE),
            ('log file',    files.DEFAULT_LOG_FILE),
        ))

        # Print such dictionary.
        stdout.output_lines(
            '{}: {}'.format(info_key, info_value)
            for info_key, info_value in info_key_to_value
        )

# --------------------( WASTELANDS                         )--------------------
# import inspect
        # Dictionary from subcommand name to _run_*() method running such
        # subcommand.
        # subcommand_name_to_method = dict(
        #     (method_name[len('_run_'):], method)
        #     for (method_name, method) in
        #         inspect.getmembers(self, inspect.ismethod)
        #     if method_name.startswith('_run_')
        # )
        # assert self._args.command_name in subcommand_name_to_method,\
        #     '"{}" not a recognized subcommand'.format(self._args.command_name)
        #
        # # Run such subcommand.
        # subcommand_name_to_method[self._args.command_name]()
        # subcommand_method = getattr(self, subcommand_method_name, None)
        # assert callable(subcommand_method),\
        #     '"{}" not callable'.format(subcommand_method_name)
        #FUXME: Implement me.
        # Run the command specified by such arguments.

        #FUXME: Display a default help message, when the user passes no
        #arguments. Does such parser already do so? This is trivial for us
        #to do as follows:
        #
        #    if len(sys.argv) == 1:
        #         self._arg_parser.print_help()
        #         return

        #FUXME: Not entirely clear as to why we require or want this dictionary.
        #Smacks of overkill, frankly. Also unclear why we require
        #_add_subparser_sim()-style methods. Just define all such subcommands
        #here, for now.

        # # Dictionary from command name to subparser object.
        # command_parsers = {}
        # command_parsers['sim'] = self._add_subparser_sim()
        #
        # # Add an identifying name and description for each command parser.
        # for command_parser_name, command_parser in command_parsers.items():
        #     command_parser.set_defaults(command_name = command_parser_name)

        # Parse command-line arguments into object attributes.
        # self._parse_args()
        # self._parse_common_args()

#FUXME; The following snippet courtesy Matthew Leingan affords an elegant means
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

# from betse.util.io import loggers
            # assert action_parser, '"{}" is None'.format(action_parser_name)
        #FUXME: This doesn't seem quite right. Don't we want to parse arguments
        #only after defining all argument syntax for actions as well?

    # @property
    # def _script_basename(self) -> str:
    #     return metadata.SCRIPT_NAME_CLI

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
