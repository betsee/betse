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

# --------------------( WASTELANDS                         )--------------------
