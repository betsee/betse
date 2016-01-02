#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
High-level application initialization common to both the CLI and GUI.
'''

# ....................{ IMPORTS                            }....................
from betse import pathtree
from betse.util.lib import libs
from betse.util.io import loggers
from betse.util.python import pythons

# ....................{ INITIALIZERS                       }....................
def init() -> None:
    '''
    Initialize the current application.

    Specifically:

    * Validate core directories and files required at program startup, creating
      all such directories and files that do *not* already exist and are
      feasible to be created.
    * Create the root logger and default handlers for such logger (e.g.,
      printing informational messages to standard output, warnings and errors to
      standard error, and appending debug messages to the default logfile).

    To support caller-specific error handling, this function is intended to be
    called immediately *after* such application begins catching otherwise
    uncaught exceptions.
    '''
    # Validate core directories and files required at program startup.
    pathtree.init()

    # Configure logging *AFTER* creating such directories, as such logging
    # writes to files in such directories.
    loggers.config.init(filename = pathtree.LOG_DEFAULT_FILENAME)
    # self._logger.error('ERROR!')
    # self._logger.warning('WARNING!')
    # self._logger.info('INFO!')
    # self._logger.debug('DEBUG!')

    # Validate mandatory dependencies *AFTER* configuring logging,
    # ensuring that exceptions raised by such validation will be logged.
    libs.init()

    # Validate the active Python interpreter *AFTER* validating mandatory
    # dependencies. While the former (mostly) comprises unenforced
    # recommendations, the latter comprises enforced requirements and hence is
    # performed first.
    pythons.init()

# --------------------( WASTELANDS                         )--------------------
    #FUXME: This appears to be required due to a PyInstaller bug. Research.
    # import tkinter.filedialog

# FUXME dear Sessums B. Delightful -- the issue with the dependancies.init() is
# still crashing the program... :(
    #FUXME: Reenable *AFTER* we integrate the CLI frontend with our
    #scientific backend. For the moment, this unacceptably provokes
    #fatal exceptions in frozen applications.
    # FUXME! This is still causing issues (AP)
