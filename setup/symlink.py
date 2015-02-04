#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''`betse`-specific `symlink` command for `setuptools`.'''

# ....................{ COMMANDS                           }....................
def add_commands(setup_options):
    '''
    Add custom `symlink` and `unsymlink` commands to the passed dictionary of
    `setuptools` options.
    '''
    #FIXME: Rethink this. A class-based approach would probably be significantly
    #more intelligible, if slightly more verbose.

    # Define such functions as closures to provide such functions read-only
    # access to such dictionary.
    def symlink():
        '''
        Editably install (e.g., in a symbolically linked manner) `betse` into
        the active Python 3 interpreter *without* performing dependency
        resolution.

        Unlike the default `develop` command, this command is suitable for
        system-wide installation.
        '''
        pass

    def unsymlink():
        '''
        Editably uninstall (e.g., in a symbolically linked manner) `betse` from
        the active Python 3 interpreter.
        '''
        pass

    setup_options['cmdclass']['symlink'] = symlink
    setup_options['cmdclass']['unsymlink'] = unsymlink

# --------------------( WASTELANDS                         )--------------------
